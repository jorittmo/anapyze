import pathlib
import os
from os.path import join,exists
import numpy as np
import pandas as pd

class ADNI(object):

    def __init__(self):

        dir_ = pathlib.Path(__file__).parent.resolve()
        self.adni_csv_dir = join(dir_, 'resources')

        amyloid_csv = join(self.adni_csv_dir, 'UCBERKELEY_AMY_6MM_13Sep2023.csv')
        csf_csv = join(self.adni_csv_dir, 'CSF.xlsx')
        asyn_csv = join(self.adni_csv_dir, 'AMPRION_ASYN_SAA_22Dec2023.csv')
        apoe_csv = join(self.adni_csv_dir, 'APOE.csv')
        adni_aging_genes_csv = join(self.adni_csv_dir, 'HS-aging_SNPs.xlsx')
        mmse_csv = join(self.adni_csv_dir, 'MMSE.csv')
        composites_csv = join(self.adni_csv_dir, 'MEM_EXEC_composites.xlsx')
        neuropsychological_csv = join(self.adni_csv_dir, 'NEUROBAT.csv')
        wmh_csv = join(self.adni_csv_dir, 'UCD_WMH.csv')

        self.amyloid_df = pd.read_csv(amyloid_csv)
        self.csf_df = pd.read_excel(csf_csv)
        self.asyn_csf_df = pd.read_csv(asyn_csv)
        self.apoe_df = pd.read_csv(apoe_csv)
        self.adni_aging_genes_df = pd.read_excel(adni_aging_genes_csv)
        self.mmse_df = pd.read_csv(mmse_csv)
        self.composites_df = pd.read_excel(composites_csv)
        self.neuropsychological_df = pd.read_csv(neuropsychological_csv)
        self.wmh_df = pd.read_csv(wmh_csv)

    @staticmethod
    def reorder_adni_data(input_dir, output_dir, dcm2niix=r'dcm2niix'):

        input_subjects = os.listdir(input_dir)

        for subject in input_subjects:
            # print(subject)
            subject_dir = join(input_dir, subject)

            for root, dirs, i_files in os.walk(subject_dir):

                for dir_ in dirs:

                    if dir_[0] == 'I':

                        input = join(root, dir_)
                        output_name = '%s_%s' % (subject, dir_)
                        output_ = join(output_dir, output_name)
                        if not exists(output_):
                            os.makedirs(output_)

                        command = dcm2niix + '-o ' + output_ + ' -f fdg ' + input
                        os.system(command)
    def is_subject_amyloid_PET_positive(self, subject_id = False, rid = False, date = 'baseline'):

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")
        df_subj = self.amyloid_df.loc[self.amyloid_df['RID'] == rid]
        print(df_subj)

        if not df_subj.empty:

            if date == 'baseline':
                ordered = df_subj.sort_values(by = ['SCANDATE'], ascending = True)
            elif date == 'last':
                ordered = df_subj.sort_values(by = ['SCANDATE'], ascending = False)
            else:
                raise Exception('Not yet implemented')

            tracer = ordered['TRACER'].values[0]
            suvr = ordered['SUMMARY_SUVR'].values[0]
            centiloids = ordered['CENTILOIDS'].values[0]
            amyloid_status = ordered['AMYLOID_STATUS'].values[0]

            return tracer, suvr, amyloid_status, centiloids

        else:
            return None, None, None, None

    def get_csf_biomarkers(self, subject_id = False, rid = False, date = 'baseline'):

        ab42 = None
        tau = None
        ptau = None
        asyn = None
        exam_date = None

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")

        df_subj = self.csf_df.loc[self.csf_df['RID'] == rid]
        df_subj['EXAMDATE'] = pd.to_datetime(df_subj.EXAMDATE)

        if not df_subj.empty:

            if date == 'baseline':
                ordered = df_subj.sort_values(by = ['EXAMDATE'], ascending = True)
            elif date == 'last':
                ordered = df_subj.sort_values(by = ['EXAMDATE'], ascending = False)
            else:
                raise Exception('Not yet implemented')

            ab42 = ordered['ABETA42'].values[0]
            tau = ordered['TAU'].values[0]
            ptau = ordered['PTAU'].values[0]

        df_subj = self.asyn_csf_df.loc[self.asyn_csf_df['RID'] == rid]

        if not df_subj.empty:

            if date == 'baseline':
                ordered = df_subj.sort_values(by = ['EXAMDATE'], ascending = True)
            elif date == 'last':
                ordered = df_subj.sort_values(by = ['EXAMDATE'], ascending = False)
            else:
                raise Exception('Not yet implemented')

            asyn_in = ordered['Result'].values[0]
            exam_date = ordered['EXAMDATE'].values[0]

            if asyn_in == 'Not_Detected':

                asyn = 0

            elif asyn_in == 'Detected-1':

                asyn = 1

            elif asyn_in == 'Detected-2':

                asyn = 2

        return ab42, tau, ptau, asyn, exam_date

    def get_genetics_data(self, subject_id = False, rid = False):

        apoe = None
        KCNMB2 = None
        TMEM106B = None
        GRN = None
        rs73069071 = None
        ABCC9 = None

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")

        filtered_apoe = self.apoe_df[self.apoe_df.RID == rid]

        if not filtered_apoe.empty:

            apoe_1 = filtered_apoe['APGEN1'].values[0]
            apoe_2 = filtered_apoe['APGEN2'].values[0]

            if apoe_1 == 4 and apoe_2 == 4:
                apoe = 2
            elif apoe_1 == 4 or apoe_2 == 4:
                apoe = 1
            else:
                apoe = 0

        genes = self.adni_aging_genes_df.loc[self.adni_aging_genes_df['RID'] == rid]

        if not genes.empty:
            KCNMB2 = genes['rs9637454_A'].values[0]
            TMEM106B = genes['rs1990622_G'].values[0]
            GRN = genes['rs5848_T'].values[0]
            rs73069071 = genes['rs73069071_C'].values[0]
            ABCC9 = genes['rs704180_G'].values[0]

        return apoe, KCNMB2, TMEM106B, GRN, rs73069071, ABCC9

    def get_cognition_data(self, subject_id = False, rid = False, date = 'baseline'):

        mmse = None
        adni_mem = None
        adni_exec = None
        adni_lan = None
        adni_vs = None

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")

        filtered_mmse = self.mmse_df[self.mmse_df.RID == rid]
        filtered_mmse['USERDATE'] = pd.to_datetime(filtered_mmse.USERDATE)

        if date == 'baseline':
            ordered = filtered_mmse.sort_values(by = ['USERDATE'], ascending = True)
        elif date == 'last':
            ordered = filtered_mmse.sort_values(by = ['USERDATE'], ascending = False)
        else:
            raise Exception('Not yet implemented')

        if not ordered.empty:
            mmse = ordered['MMSCORE'].values[0]

        filtered_composite = self.composites_df[self.composites_df.RID == rid]
        filtered_composite['EXAMDATE'] = pd.to_datetime(filtered_composite.EXAMDATE)

        if date == 'baseline':
            ordered = filtered_composite.sort_values(by = ['EXAMDATE'], ascending = True)
        elif date == 'last':
            ordered = filtered_composite.sort_values(by = ['EXAMDATE'], ascending = False)
        else:
            raise Exception('Not yet implemented')

        if not ordered.empty:

            if not ordered['ADNI_MEM'].dropna().empty:
                adni_mem = ordered['ADNI_MEM'].dropna().values[0]

            if not ordered['ADNI_EF'].dropna().empty:
                adni_exec = ordered['ADNI_EF'].dropna().values[0]

            if not ordered['ADNI_LAN'].dropna().empty:
                adni_lan = ordered['ADNI_LAN'].dropna().values[0]

            if not ordered['ADNI_VS'].dropna().empty:
                adni_vs = ordered['ADNI_VS'].dropna().values[0]

        return mmse, adni_mem, adni_exec, adni_lan, adni_vs

    def get_mmse_closest_to_date(self, date, subject_id = False, rid = False):

        from datetime import datetime
        from dateutil.relativedelta import relativedelta

        #convert string to datetime
        date = datetime.strptime(date, '%Y-%m-%d')

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")

        filtered_mmse = self.mmse_df[self.mmse_df.RID == rid]
        filtered_mmse["USERDATE"] = pd.to_datetime(filtered_mmse.USERDATE)

        diff_time_mmse = 1000
        mmse = 30
        date_out = '1900-01-01'

        for mmse_index, mmse_row in filtered_mmse.iterrows():

            date_mmse = mmse_row["USERDATE"]
            mmse_ = mmse_row["MMSCORE"]

            diff_time_ = (
                    relativedelta(date_mmse, date).years
                    + relativedelta(date_mmse, date).months / 12
                    + relativedelta(date_mmse, date).days / 365
            )

            if np.abs(diff_time_) < diff_time_mmse:
                diff_time_mmse = np.abs(diff_time_)
                mmse = mmse_
                date_out = date_mmse

        return mmse, date_out

    def get_neuropsychological_battery(self, subject_id = False, rid = False, date = 'baseline'):

        # REF: https://adni.bitbucket.io/reference/neurobat.html

        clock_drawing = None
        clock_copy = None
        fluency_animals = None
        tmt_a = None
        tmt_b = None
        boston_naming_spontaneus = None
        boston_naming_total = None
        delayed_recall = None

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")

        filtered_npb = self.neuropsychological_df[self.neuropsychological_df.RID == rid]
        filtered_npb['USERDATE'] = pd.to_datetime(filtered_npb.USERDATE)

        if date == 'baseline':
            ordered = filtered_npb.sort_values(by = ['USERDATE'], ascending = True)
        elif date == 'last':
            ordered = filtered_npb.sort_values(by = ['USERDATE'], ascending = False)
        else:
            raise Exception('Not yet implemented')

        if not ordered.empty:

            if not ordered['CLOCKSCOR'].dropna().empty:
                clock_drawing = ordered['CLOCKSCOR'].dropna().values[0]

            if not ordered['COPYSCOR'].dropna().empty:
                clock_copy = ordered['COPYSCOR'].dropna().values[0]

            if not ordered['CATANIMSC'].dropna().empty:
                fluency_animals = ordered['CATANIMSC'].dropna().values[0]

            if not ordered['TRAASCOR'].dropna().empty:
                tmt_a = ordered['TRAASCOR'].dropna().values[0]

            if not ordered['TRABSCOR'].dropna().empty:
                tmt_b = ordered['TRABSCOR'].dropna().values[0]

            if not ordered['BNTSPONT'].dropna().empty:
                boston_naming_spontaneus = ordered['BNTSPONT'].dropna().values[0]

            if not ordered['BNTTOTAL'].dropna().empty:
                boston_naming_total = ordered['BNTTOTAL'].dropna().values[0]

            if not ordered['AVDEL30MIN'].dropna().empty:
                delayed_recall = ordered['AVDEL30MIN'].dropna().values[0]

        return (clock_drawing, clock_copy, fluency_animals, tmt_a, tmt_b,
                boston_naming_spontaneus, boston_naming_total, delayed_recall)

    def get_wmh(self, subject_id = False, rid = False, date = 'baseline'):

        wmh_vol = None
        left_hippo_vol = None
        right_hippo_vol = None
        tiv = None

        if subject_id:
            rid = int(subject_id[-4:])
        elif not subject_id and not rid:
            raise ValueError("Subject_ID or RID must be provided")

        filtered_wmh = self.wmh_df[self.wmh_df.RID == rid]
        filtered_wmh['EXAMDATE'] = pd.to_datetime(filtered_wmh.EXAMDATE)

        if date == 'baseline':
            ordered = filtered_wmh.sort_values(by = ['EXAMDATE'], ascending = True)
        elif date == 'last':
            ordered = filtered_wmh.sort_values(by = ['EXAMDATE'], ascending = False)
        else:
            raise Exception('Not yet implemented')

        if not ordered.empty:

            if not ordered['TOTAL_WMH'].dropna().empty:
                wmh_vol = ordered['TOTAL_WMH'].dropna().values[0]

            if not ordered['LEFT_HIPPO'].dropna().empty:
                left_hippo_vol = ordered['LEFT_HIPPO'].dropna().values[0]

            if not ordered['RIGHT_HIPPO'].dropna().empty:
                right_hippo_vol = ordered['RIGHT_HIPPO'].dropna().values[0]

            if not ordered['TOTAL_BRAIN'].dropna().empty:
                tiv = ordered['TOTAL_BRAIN'].dropna().values[0]

        return wmh_vol, left_hippo_vol, right_hippo_vol, tiv