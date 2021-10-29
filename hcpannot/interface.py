################################################################################
# interface.py
#
# Interface for drawing the HCP lines.
# by Noah C. Benson <nben@uw.edu>

# Import things
import sys, os, pimms
import numpy as np
import pyrsistent as pyr
import matplotlib as mpl
import matplotlib.pyplot as plt
import ipywidgets as widgets
import neuropythy as ny

from .core import image_order

# The default image size; this is the assumed size of the images that are
# downloaded from the OSF and displayed, in pixels.
default_imshape = (864*2, 864*2)
# The default x- and y-value limits of the flatmaps that were used for
# generating the images.
default_xlim = (-100, 100)
default_ylim = (-100, 100)
# The default grid used for the display. The None is a stand-in for the image of
# the current contour's highlight.
default_grid = ((None,        'polar_angle'),
                ('curvature', 'eccentricity'))
# The path that we load images from by default.
default_load_path = '/data'
default_osf_url = 'osf://tery8/'

def imgrid_to_xy(pts,
                 grid=default_grid,
                 imshape=default_imshape,
                 xlim=default_xlim,
                 ylim=default_ylim):
    '''
    `imgrid_to_xy(pts)` yields a 2xN matrix the same size as the given
      (2xN) matrix `pts`, for which the points have been converted from
      coordinates in the given image grid (`grid` option).
    '''
    pts = np.array(pts)
    (R,C) = imshape[:2]
    rg = len(grid)
    cg = len(grid[0])
    (rs,cs) = (R/rg, C/cg)
    rmu = (rs-1)/2
    cmu = (cs-1)/2
    (xmin,xmax) = xlim
    (ymin,ymax) = ylim
    xmu = 0.5*(xmin + xmax)
    ymu = 0.5*(ymin + ymax)
    rpx2yu = -(ymax - ymin) / r
    cpx2xu = (xmax - xmin) / c
    (c,r) = pts if pts.shape[0] == 2 else pts.T
    while True:
        ii = c > cs
        if len(ii) == 0: break
        c[ii] -= cs
    while True:
        ii = r > rs
        if len(ii) == 0: break
        r[ii] -= rs
    x = xmu + (cs - cmu)*cpx2xu
    y = ymu + (rs - rmu)*rpx2yu
    return np.array([x,y])
def xy_to_imgrid(pts,
                 grid=default_grid,
                 imshape=default_imshape,
                 xlim=default_xlim,
                 ylim=default_ylim):
    '''
    `xy_to_imgrid(pts)` yields a 2xN matrix the same size as the given
      (2xN) matrix `pts`, for which the points have been converted from
      coordinates in the default flatmap representation to the given
      image grid (`grid` option).
    '''
    pts = np.asarray(pts)
    (R,C) = imshape[:2]
    rg = len(grid)
    cg = len(grid[0])
    (rs,cs) = (R/rg, C/cg)
    rmu = (rs-1)/2
    cmu = (cs-1)/2
    (xmin,xmax) = xlim
    (ymin,ymax) = ylim
    xmu = 0.5*(xmin + xmax)
    ymu = 0.5*(ymin + ymax)
    yu2rpx = -rs / (ymax - ymin)
    xu2cpx = cs / (xmax - xmin)
    (x,y) = pts if pts.shape[0] == 2 else pts.T
    c = cmu + (x - xmu)*xu2cpx
    r = rmu + (y - ymu)*yu2rpx
    return np.array([[[c+c0*cs,r+r0*rs] for (c0,_) in enumerate(g)]
                     for (r0,g) in enumerate(grid)])

def point_decorate_plot(ax, pts, *args, **kw):
    grid = kw.pop('grid', default_grid)
    imshape = kw.pop('imshape', default_imshape)
    xlim = kw.pop('xlim', default_xlim)
    ylim = kw.pop('ylim', default_ylim)
    rcs = xy_to_imgrid(pts, grid=grid, imshape=imshape, xlim=xlim, ylim=ylim)
    plots = [ax.plot(c, r, *args, **kw)
             for row in rcs
             for (r,c) in row]
    return plots
def segs_decorate_plot(ax, segs, *args, **kw):
    from matplotlib.collections import LineCollection as lncol
    grid = kw.pop('grid', default_grid)
    imshape = kw.pop('imshape', default_imshape)
    xlim = kw.pop('xlim', default_xlim)
    ylim = kw.pop('ylim', default_ylim)
    pts_nx2 = np.reshape(segs, (-1, 2))
    rcs = xy_to_imgrid(pts_nx2.T, grid=grid, imshape=imshape, xlim=xlim, ylim=ylim)
    plots = [lncol(segs, *args, **kw)
             for row in rcs
             for segs0 in row
             for segs in [np.reshape(segs0.T, (-1,2,2))]]
    for p in plots:
        ax.add_collection(p)
    return plots
def clicks_decorate_plot(ax, pts, *args, **kw):
    grid = kw.pop('grid', default_grid)
    imshape = kw.pop('imshape', default_imshape)
    xlim = kw.pop('xlim', default_xlim)
    ylim = kw.pop('ylim', default_ylim)
    (rs,cs) = imshape[:2]
    rs /= len(grid)
    cs /= len(grid[0])
    if len(pts) > 0:
        (x,y) = np.transpose(pts)
    else:
        (x,y) = ([], [])
    x = np.mod(x, cs)
    y = np.mod(y, rs)
    plots = []
    for (r,row) in enumerate(grid):
        for (c,col) in enumerate(row):
            pp = ax.plot(x + cs*c, y + rs*r, *args, **kw)
            for p in pp: plots.append(p)
    return plots
def clicks_update_plot(ax, plots, pts, grid=default_grid, imshape=default_imshape):
    (rs,cs) = imshape[:2]
    rs /= len(grid)
    cs /= len(grid[0])
    (x,y) = np.transpose(pts)
    x = np.mod(x, cs)
    y = np.mod(y, rs)
    for plot in plots:
        (px, py) = plot.get_data()
        if len(px) > 0:
            dx = px[0] - x[0]
            dy = py[0] - y[0]
            plot.set_data(x+dx, y+dy)
    return plots
def load_subimage(sid, h, name,
                  load_path=default_load_path, osf_url=default_osf_url):
    from PIL import Image
    flnm = os.path.join(load_path, str(sid), '%d_%s_%s.png' % (sid, h, name))
    with Image.open(flnm) as im:
        arr = np.array(im)
    return arr
def curry_load_subimage(sid, h, name,
                        load_path=default_load_path, osf_url=default_osf_url):
    return lambda:load_subimage(sid, h, name,
                                load_path=load_path, osf_url=osf_url)
def load_subwang(sid, h, load_path=default_load_path, osf_url=default_osf_url):
    import neuropythy as ny
    flnm = os.path.join(load_path, str(sid), '%d_%s_wang.mgz' % (sid, h))
    return np.array(ny.load(flnm, 'mgh', to='data'))
def imcat(grid):
    col = [np.concatenate(row, axis=1) for row in grid]
    return np.concatenate(col, axis=0)
def plot_imcat(ims, grid, k):
    grid = [[ims[k if g is None else g] for g in row]
            for row in grid]
    return imcat(grid)
def prep_subdata(sid, h, load_path=default_load_path, osf_url=default_osf_url):
    dirname = os.path.join(load_path, str(sid))
    if not os.path.isfile(dirname):
        pp = ny.util.pseudo_path(osf_url)
        path = pp.local_path('annot-images', '%d.tar.gz' % sid)
        import tarfile
        with tarfile.open(path) as fl:
            fl.extractall(load_path)
    ims = {imname: curry_load_subimage(sid, h, imname,
                                       load_path=load_path, osf_url=osf_url)
           for imname in image_order}
    ims['wang'] = lambda:load_subwang(sid, h,
                                      load_path=load_path, osf_url=osf_url)
    return pimms.lmap(ims)
def curry_prep_subdata(sid, h,
                       load_path=default_load_path, osf_url=default_osf_url):
    return lambda:prep_subdata(sid, h, load_path=load_path, osf_url=osf_url)

subject_ids = (100610, 102311, 102816, 104416, 105923, 108323, 109123, 111312,
               111514, 114823, 115017, 115825, 116726, 118225, 125525, 126426,
               128935, 130114, 130518, 131217, 131722, 132118, 134627, 134829,
               135124, 137128, 140117, 144226, 145834, 146129, 146432, 146735,
               146937, 148133, 150423, 155938, 156334, 157336, 158035, 158136,
               159239, 162935, 164131, 164636, 165436, 167036, 167440, 169040,
               169343, 169444, 169747, 171633, 172130, 173334, 175237, 176542,
               177140, 177645, 177746, 178142, 178243, 178647, 180533, 181232,
               181636, 182436, 182739, 185442, 186949, 187345, 191033, 191336,
               191841, 192439, 192641, 193845, 195041, 196144, 197348, 198653,
               199655, 200210, 200311, 200614, 201515, 203418, 204521, 205220,
               209228, 212419, 214019, 214524, 221319, 233326, 239136, 246133,
               249947, 251833, 257845, 263436, 283543, 318637, 320826, 330324,
               346137, 352738, 360030, 365343, 380036, 381038, 385046, 389357,
               393247, 395756, 397760, 401422, 406836, 412528, 429040, 436845,
               463040, 467351, 525541, 536647, 541943, 547046, 550439, 552241,
               562345, 572045, 573249, 581450, 585256, 601127, 617748, 627549,
               638049, 644246, 654552, 671855, 680957, 690152, 706040, 724446,
               725751, 732243, 751550, 757764, 765864, 770352, 771354, 782561,
               783462, 789373, 814649, 818859, 825048, 826353, 833249, 859671,
               861456, 871762, 872764, 878776, 878877, 898176, 899885, 901139,
               901442, 905147, 910241, 926862, 927359, 942658, 943862, 951457,
               958976, 966975, 971160, 973770, 995174)

subject_data = pimms.lmap({(sid,h): curry_prep_subdata(sid, h)
                           for sid in subject_ids
                           for h in ['lh','rh']})

boundary_contours = pyr.pmap(
    {'V1-middle': 'isoang_90',
     'V1/V2 dorsal': 'isoang_vml',
     'V1/V2 ventral': 'isoang_vmu',
     'V2/V3 dorsal': 'isoang_hml',
     'V2/V3 ventral': 'isoang_hmu',
     'V3/Outer dorsal': 'isoang_vml',
     'V3/Outer ventral': 'isoang_vmu'})
contour_names = tuple(['0.5° iso-eccen'] + 
                      ['%d° iso-eccen' % k for k in [1,2,4,7]] +
                      list(boundary_contours.keys()))
contour_key = dict(boundary_contours)
contour_key['0.5° iso-eccen'] = 'isoecc_0.5'
for k in [1,2,4,7]:
    contour_key['%d° iso-eccen' % k] = 'isoecc_%d' % k
contour_key = pyr.pmap(contour_key)

# #ROITool #####################################################################
class ROITool(object):
    '''
    ROITool is a tool for drawing ROIs and contours on the HCP.
    '''
    def __init__(self,
                 figsize=2.5, sidepanel_width='250px', dropdown_width='85%',
                 start_contour='V3/Outer ventral',
                 imshape=default_imshape, grid=default_grid):
        self.imshape = imshape
        self.grid = grid
        self.start_contour = start_contour
        start_contour = contour_key[start_contour]
        (grid_rs, grid_cs) = (len(grid), len(grid[0]))
        figh = figsize * grid_rs / grid_cs
        # We store clicks in auto-dicts.
        self.clicks = {sid: {h: ny.util.auto_dict(None, [])
                             for h in ['lh','rh']}
                       for sid in subject_ids}
        # Go ahead and setup the Widgets.
        self.sid_select = widgets.Dropdown(
            options=subject_ids,
            value=subject_ids[0],
            description='SID:',
            layout={'width': dropdown_width})
        self.hemi_select = widgets.Dropdown(
            options=['LH','RH'],
            value='LH',
            description='Hemi:',
            layout={'width': dropdown_width})
        self.line_select = widgets.Dropdown(
            options=contour_names,
            value=self.start_contour,
            description='Contour:',
            layout={'width': dropdown_width})
        self.anat_shown = widgets.Checkbox(
            description='Wang et al. (2015) Contours',
            value=True)
        self.anat_color = widgets.ColorPicker(
            description='Wang Color:',
            concise=True,
            value='white',
            layout={'width':'50%'})
        self.contour_shown = widgets.Checkbox(
            description='Drawn Contour',
            value=True)
        self.contour_color = widgets.ColorPicker(
            description='Draw Color:',
            concise=True,
            value='cyan',
            layout={'width':'50%'})
        self.notes_area = widgets.Textarea(
            value='', 
            description='',
            layout={'height': '%dpx' % (200), 'width': '95%'})
        self.notes_panel = widgets.VBox(
            [widgets.Label('Contour Notes:'), self.notes_area],
            layout={'align_items': 'flex-start', 'width':'100%', 'height':'200px'})
        self.save_button = widgets.Button(description='Save.')
        #center_layout = widgets.Layout(align_items='center')
        self.save_box = widgets.HBox(
            children=[self.save_button],
            layout={'align_items': 'center'})
        self.controls = (self.sid_select,
                         self.hemi_select,
                         self.line_select,
                         self.anat_shown,
                         self.anat_color,
                         self.contour_shown,
                         self.contour_color,
                         self.notes_panel,
                         self.save_button)
        self.control_panel = widgets.VBox(self.controls, layout={'height': '95%'})
        # The start/default values:
        self.sid = self.sid_select.value
        self.hemi = self.hemi_select.value.lower()
        # Setup the figure.
        subdata = subject_data[(self.sid, self.hemi)]
        im0 = plot_imcat(subdata, grid, start_contour)
        segs = subdata['wang']
        #(im_rs, im_cs) = imshape
        (dot_rs, dot_cs) = imshape #(im_rs*grid_rs, im_cs*grid_cs)
        (fig,ax) = plt.subplots(
            constrained_layout=True,
            figsize=(figsize, figsize*grid_rs/grid_cs),
            dpi=72*4)
        self.figure = fig
        self.axes = ax
        fig.canvas.toolbar_position = 'bottom'
        fig.canvas.title_visible = False
        fig.canvas.header_visible = False
        fig.canvas.footer_visible = False
        #ax.format_coord = lambda x,y: ''
        self.wang_plot = segs_decorate_plot(
            ax, segs, color=self.anat_color.value, lw=0.3, zorder=10,
            grid=grid, imshape=imshape)
        self.contour_plot = []
        # Initialize the display for this subject/hemi
        self.image_plot = ax.imshow(im0)
        ax.axis('off')
        # Setup the listener functions...
        self.sid_select.observe(ny.util.curry(self.update, 'sid'), 'value')
        self.hemi_select.observe(ny.util.curry(self.update, 'hemi'), 'value')
        self.line_select.observe(ny.util.curry(self.update, 'line'), 'value')
        self.anat_shown.observe(ny.util.curry(self.update, 'anat'), 'value')
        self.contour_shown.observe(ny.util.curry(self.update, 'contour'), 'value')
        self.anat_color.observe(ny.util.curry(self.update, 'anatcolor'), 'value')
        self.contour_color.observe(ny.util.curry(self.update, 'contourcolor'), 'value')
        self.canvas_conns = [
            #fig.canvas.mpl_connect('close_event', self.on_close),
            fig.canvas.mpl_connect('button_press_event', self.on_click)]
        # Final touches:
        self.control_panel.layout = widgets.Layout(width=sidepanel_width,
                                                   height='600px',
                                                   align_items='center')
        pane = widgets.HBox(
            [self.control_panel, fig.canvas],
            layout=widgets.Layout(
                flex_flow='row',
                align_items='center',
                width='100%',
                height='700px',
                border='#000000'))
        display(pane)
        # For saving errors that get caught in events:
        self._event_error = None

    def update(self, var, change):
        sid = int(self.sid_select.value)
        h = self.hemi_select.value.lower()
        contour = self.line_select.value
        anat = self.anat_shown.value
        ax = self.axes
        fig = self.figure
        # Get & remove the (now deprecated) plots:
        implot = self.image_plot
        wangplot = self.wang_plot
        # What updated?
        if var == 'sid':
            # What's the new control selection:
            sid = int(change.new)
            h = self.hemi
            # Remove contour plots if need-be
            for c in self.contour_plot:
                c.remove()
            self.contour_plot = []
            # New plots:
            subdata = subject_data[(sid, h)]
            contour = self.start_contour
            c = contour_key[contour]
            im0 = plot_imcat(subdata, grid, c)
            self.image_plot.set_data(im0)
            for ln in self.wang_plot: ln.remove()
            segs = subdata['wang']
            self.wang_plot = segs_decorate_plot(
                ax, segs,
                grid=self.grid, imshape=self.imshape,
                color=self.anat_color.value, lw=0.3, zorder=10)
            pts = self.clicks[sid][h][contour]
            self.contour_plot = clicks_decorate_plot(
                    ax, pts, 'o-',
                    grid=self.grid, imshape=self.imshape,
                    color=self.contour_color.value, lw=0.5, ms=1.5)
            # Update the output data:
            self.sid = sid
            # Update the controls:
            #anat_shown.value = True
            self.line_select.value = self.start_contour
            self.contour_shown.value = True
        elif var == 'hemi':
            # New plots:
            h = change.new
            subdata = subject_data[(self.sid, h)]
            contour = self.start_contour
            c = contour_key[contour]
            im0 = plot_imcat(subdata, grid, c)
            self.image_plot.set_data(im0)
            # Update Wang plot lines:
            for ln in self.wang_plot: ln.remove()
            segs = subdata['wang']
            self.wang_plot = segs_decorate_plot(
                ax, segs,
                grid=self.grid, imshape=self.imshape,
                color=self.anat_color.value, lw=0.3, zorder=10)
            for ln in self.wang_plot:
                ln.set_visible(self.anat_shown.value)
            # And the drawn contours:
            for c in self.contour_plot: c.remove()
            pts = self.clicks[sid][h][contour]
            self.contour_plot = clicks_decorate_plot(
                    ax, pts, 'o-',
                    grid=self.grid, imshape=self.imshape,
                    color=self.contour_color.value, lw=0.5, ms=1.5)
            # Update the output data:
            self.hemi = h
            # Update the controls:
            #anat_shown.value = True
            self.line_select.value = self.start_contour
            self.contour_shown.value = True
        elif var == 'line':
            # Remove contour plots if need-be
            for c in self.contour_plot: c.remove()
            contour = change.new
            c = contour_key[contour]
            subdata = subject_data[(sid,h)]
            im = plot_imcat(subdata, self.grid, c)
            self.image_plot.set_data(im)
            self.contour_shown.value = True
            # Redraw the chosen contours if need-be
            pts = self.clicks[sid][h][contour]
            self.contour_plot = clicks_decorate_plot(
                ax, pts, 'o-',
                grid=self.grid, imshape=self.imshape,
                color=self.contour_color.value, lw=0.5, ms=1.5)
            for ln in self.contour_plot:
                ln.set_visible(True)
        elif var == 'anat':
            anat = change.new
            for ln in self.wang_plot: ln.set_visible(anat)
        elif var == 'contour':
            c = change.new
            cplot = self.contour_plot
            for ln in cplot: ln.set_visible(c)
        elif var == 'anatcolor':
            c = change.new
            for ln in self.wang_plot: ln.set_color(c)
        elif var == 'contourcolor':
            c = change.new
            cplot = self.contour_plot
            for ln in cplot: ln.set_color(c)
        else: return
        fig.canvas.draw_idle()
        return None
    
    # Setup the figure clicks!
    def on_click(self, event):
        try:
            ax = self.axes
            fig = self.figure
            if event.inaxes != ax: return
            sid = int(self.sid_select.value)
            h = self.hemi_select.value.lower()
            contour = self.line_select.value
            c = contour_key[contour]
            cplot = self.contour_plot
            pts = self.clicks[sid][h][contour]
            # if shift is down, we delete the last point
            if event.key == 'shift':
                if len(pts) == 0: return
                del pts[-1]
                if len(pts) == 0:
                    self.contour_plot = []
                    for ln in cplot: ln.remove()
                    fig.canvas.draw()
                    return
            else: # add the points
                pts.append((event.xdata, event.ydata))
            # redraw the line regardless
            if len(pts) == 1:
                self.contour_plot = clicks_decorate_plot(
                    ax, pts, 'o-',
                    grid=self.grid, imshape=self.imshape,
                    color=self.contour_color.value, lw=0.5, ms=1.5)
                for ln in self.contour_plot:
                    ln.set_visible(True)
            else:
                clicks_update_plot(
                    ax, self.contour_plot, pts,
                    grid=self.grid, imshape=self.imshape)
            fig.canvas.draw()
        except Exception as e:
            self._event_error = sys.exc_info()
            raise
        
        
        
