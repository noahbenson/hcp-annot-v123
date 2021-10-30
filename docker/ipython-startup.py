

import hcpannot
from hcpannot import ROITool

# A hack to pull the github username out of the git directory
# and put it where ROITool can find it.
import os
git_username = os.popen("grep 'url = https://github.com/' ~/.extgit/config").read()
git_username = git_username.split('github.com/')[-1].split('/hcp-annot-v123')[0]
os.environ['GIT_USERNAME'] = git_username





    
