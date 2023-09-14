import pathlib
import pandas as pd
from os.path import join


class adni(object):

    def __init__(self):

        dir_ = pathlib.Path(__file__).parent.resolve()
        self.adni_csv_dir = join(dir_,'resources', 'adni')

        amyloid_csv = join(self.adni_csv_dir, 'UCBERKELEY_AMY_6MM_13Sep2023.csv')

        csf_csv = join(self.adni_csv_dir, 'CSF.xlsx')


        self.amyloid_df = pd.read_csv(amyloid_csv)
        self.csf_df = pd.read_excel(csf_csv)

    def is_subject_amyloid_positive(self, subject_id=False, rid=False, date='baseline'):

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")
        df_subj = self.amyloid_df.loc[self.amyloid_df['RID'] == rid]
        print(df_subj)

        if not df_subj.empty:

            if date == 'baseline':
                ordered = df_subj.sort_values(by=['SCANDATE'], ascending=True)
            elif date == 'last':
                ordered = df_subj.sort_values(by=['SCANDATE'], ascending=False)
            else:
                raise Exception('Not yet implemented')

            tracer = ordered['TRACER'].values[0]
            suvr = ordered['SUMMARY_SUVR'].values[0]
            centiloids = ordered['CENTILOIDS'].values[0]
            amyloid_status = ordered['AMYLOID_STATUS'].values[0]

            return tracer, suvr, amyloid_status, centiloids

        else:
            return False, False, False, False

    def get_csf_biomarkers(self,subject_id=False, rid=False, date='baseline'):

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")

        df_subj = self.csf_df.loc[self.csf_df['RID'] == rid]

        if not df_subj.empty:

            if date == 'baseline':
                ordered = df_subj.sort_values(by=['EXAMDATE'], ascending=True)
            elif date == 'last':
                ordered = df_subj.sort_values(by=['EXAMDATE'], ascending=False)
            else:
                raise Exception('Not yet implemented')

            ab42 = ordered['ABETA42'].values[0]
            tau = ordered['TAU'].values[0]
            ptau = ordered['PTAU'].values[0]

            print(ab42, tau, ptau)

            return ab42, tau, ptau

        else:
            return False, False, False


