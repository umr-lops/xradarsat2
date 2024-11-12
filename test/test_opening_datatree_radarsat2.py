import pdb
from xradarsat2.utils import get_glob,load_config
import xradarsat2
import time
import logging
logging.basicConfig(level=logging.DEBUG)
logging.debug('start opening RadarSAT-2 product')
# conf = getconfig.get_config()
# subswath = conf['product_paths'][0]

t0 = time.time()
conf = load_config()
folder_path = conf['folder_path']
dt = xradarsat2.rs2_reader(folder_path)
elapse_t = time.time()-t0

print(type(dt),dt)
print('out of the reader')
print(dt)
print('time to read the SAFE through nfs: %1.2f sec'%elapse_t)
dt = xradarsat2.load_digital_number(dt,chunks={'pol':'VV','line':6000,'sample':8000})
print('DN',dt['/digital_numbers'])
# pdb.set_trace()
