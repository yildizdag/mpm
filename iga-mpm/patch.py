import numpy as np
from readPatch2D import *
from nurbs import *

class patch(readPatch2D,nurbs):
    def __init__(self,file_ID):
        super().__init__(file_ID)
