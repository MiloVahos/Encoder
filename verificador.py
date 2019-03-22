# -*- coding: utf-8 -*-
import filecmp

print filecmp.cmp('./Preambulos.txt','./PreambulosEst.txt',shallow=False)