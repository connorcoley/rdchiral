try:
	import rdkit
except:
	raise ImportError('RDKit cannot be imported. rdchiral is a wrapper for RDKit, instructions for installing RDKit can be found at https://www.rdkit.org/docs/Install.html')
if rdkit.__version__ < '2019':
	raise ImportWarning('RDKit version {} may not be compatable with rdchiral. Upgrade to RDKit version 2019 or higher.'.format(rdkit.__version__))
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL) 