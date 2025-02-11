import adni_tools

dir_data_ = '/home/jsilva/data/ADNI_data'

proc = adni_tools.ADNI()

proc.reorder_adni_data(dir_data_)