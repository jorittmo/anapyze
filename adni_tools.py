import pathlib
import pandas as pd
from os.path import join, exists


class adni(object):

    def __init__(self):

        dir_ = pathlib.Path(__file__).parent.resolve()
        self.adni_csv_dir = join(dir_,'resources', 'ADNI_CSVs')

        amyloid_csv = join(self.adni_csv_dir, 'UCBERKELEY_AMY_6MM_12Jul2023.csv')
        print('Amyloid csv = %s' % amyloid_csv)

        self.amyloid_df = pd.read_csv(amyloid_csv)

    def is_subject_amyloid_positive(self, subject_id, date='baseline'):

        rid = int(subject_id[-4:])
        df_subj = self.amyloid_df.loc[self.amyloid_df['RID'] == rid]

        if not df_subj.empty:

            if date == 'baseline':
                ordered = df_subj.sort_values(by=['EXAMDATE'], ascending=True)
            elif date == 'last':
                ordered = df_subj.sort_values(by=['EXAMDATE'], ascending=False)
            else:
                raise Exception('Not yet implemented')

            tracer = ordered['TRACER'].values[0]
            suvr = ordered['SUMMARY_SUVR'].values[0]
            centiloids = ordered['CENTILOIDS'].values[0]
            amyloid_status = ordered['AMYLOID_STATUS'].values[0]

            return tracer, suvr, amyloid_status, centiloids

        else:
            return False
