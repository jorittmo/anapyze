import pathlib
import os
from os.path import join, exists
import numpy as np
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta


class ADNI(object):
    """
    This class provides utilities for handling and processing ADNI data.
    """

    def __init__(self):
        """
        Initializes the ADNI class with the necessary directories and file paths.
        """

        dir_repos = pathlib.Path(__file__).parents[1]
        print(dir_repos)

        self.adni_csv_dir = join(dir_repos, 'adni_resources', 'resources')

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
    def reorder_adni_data(input_dir, output_dir, dcm2niix = r'dcm2niix'):
        """
        Reorders the ADNI data based on the input directory and output directory.

        :param input_dir: The directory where the input data is located.
        :param output_dir: The directory where the reordered data will be saved.
        :param dcm2niix: The path to the dcm2niix tool.
        """

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

    @staticmethod
    def filter_mri_csv(mri_csv, output_csv):
        """
        Filters the MRI CSV file and saves the result in a new CSV file.

        :param mri_csv: The path to the MRI CSV file.
        :param output_csv: The path to the output CSV file.
        """
        mri_df = pd.read_csv(mri_csv)

        filtered_mri_df = pd.DataFrame(
                columns = ['SUBJECT_ID', 'MRI_AGE', 'MRI_ID', 'DESCRIPTION']
                )

        out = 'test.csv'

        subjects = mri_df['Subject ID'].unique()
        for subject in subjects:

            subject_df = mri_df.loc[mri_df['Subject ID'] == subject]

            ages = subject_df['Age'].unique()

            for age in ages:

                print(subject, age)

                age_df = subject_df.loc[subject_df['Age'] == age]

                age_df = age_df.loc[age_df['Description'] != '3 Plane Localizer']
                age_df = age_df.loc[age_df['Description'] != '3_Plane_Localizer']
                age_df = age_df.loc[age_df['Description'] != 'localizer']
                age_df = age_df.loc[age_df['Description'] != 'Calibration Scan']
                age_df = age_df.loc[age_df['Description'] != 'Field Mapping']
                age_df = age_df.loc[age_df['Description'] != 'Axial Field Mapping']
                age_df = age_df.loc[age_df['Description'] != 'B1-Calibration PA']
                age_df = age_df.loc[age_df['Description'] != 'B1-Calibration Body']

                if age_df.empty:
                    print('No MRI for subject %s' % subject)
                    continue
                else:
                    if len(age_df) == 1:
                        id = age_df['Image ID'].values[0]
                        description = age_df['Description'].values[0]
                        new_row = pd.DataFrame(
                                [{'SUBJECT_ID': subject, 'MRI_AGE': age, 'MRI_ID': id, 'DESCRIPTION': description}],
                                columns = ['SUBJECT_ID', 'MRI_AGE', 'MRI_ID', 'DESCRIPTION']
                                )
                        filtered_mri_df = pd.concat([filtered_mri_df, new_row], ignore_index = True)

                    else:
                        descriptions = []
                        for mri_index, mri_row in age_df.iterrows():
                            this_image_id = mri_row['Image ID']
                            this_description = mri_row['Description']

                            print('%s: %s' % (this_image_id, this_description))

                            descriptions.append(this_description)

                        accepted_descriptions = set(
                                ['Accelerated SAG IR-SPGR', 'MPRAGE SENSE2', 'MPRAGE GRAPPA2',
                                 'MPRAGE_GRAPPA2', 'Accelerated Sag IR-FSPGR', 'IR-SPGR w/acceleration'
                                                                               'Accelerated Sagittal MPRAGE',
                                 'Accelerated SAG IR-FSPGR',
                                 'Accelerated Sag IR-SPGR', 'Accelerated Sagittal MPRAGE',
                                 'MPRAGE SENSE2 SENSE', 'MPRAGE SENSE', 'Sagittal 3D Accelerated MPRAGE',
                                 'Sagittal 3D Accelerated 0 angle MPRAGE', 'ACCELERATED SAG IR-SPGR',
                                 'Accelerated Sagittal IR-FSPGR', 'IR-SPGR w/acceleration',
                                 'Accelerated Sagittal MPRAGE Phase A-P', 'Sag Accel IR-FSPGR',
                                 'Accelerated_Sag_IR-FSPGR', 'MPRAGE_GRAPPA2          straight no angle',
                                 'MPRAGE SENS', 'MPR; GradWarp; B1 Correction; N3; Scaled <- MP-RAGE',
                                 'MPRAGE REPEAT', 'SAG IR-FSPGR-Repeat', 'IR-FSPGR-Repeat',
                                 'IR-FSPGR REPEAT', 'MPRAGE_ASO_repeat', 'MPRAGE_ Sag  - NO ANGLE=',
                                 'MPRAGE Repeat', 'MP-RAGE REPEAT', 'MP RAGE SAGITTAL REPEAT',
                                 'MP-RAGE-REPEAT', 'MP-RAGE-Repeat', 'MP-RAGE repeat',
                                 'MP-RAGE', 'Sag IR-FSPGR Repeat', 'MP RAGE REPEAT',
                                 'ADNI       MPRAGE', 'MPRAGE', 'ASO-MPRAGE', 'MPRAGE AUTOSHIM ON',
                                 'Sag IR-SPGR', 'MPRAGE SAGITTAL', 'MPRAGE ASO', 'ADNI SH    MPRAGE ASO',
                                 'MPRAGE SAG', 'ADNI-R11   MPRAGE-REPEA', 'ADNI-R11   MPRAGE',
                                 'ADNI-R11-ASASO-MPRAGE', 'ADNI-R11-ASASO-MPRAGE(2',
                                 'ADNI       MPRAGEASOREP', 'ADNI       MPRAGE ASO', 'MPRAGEREPEATASO',
                                 'MPRAGEASO', 'REPEAT SAG 3D MPRAGE', 'SAG 3D MPRAGE',
                                 'REPEAT SAG 3D MP RAGE NO ANGLE', 'SAG 3D MPRAGE NO ANGLE',
                                 'SAG MPRAGE GRAPPA2 NO ANGLE', 'SAG MPRAGE NO ANGLE',
                                 'SAG MP-RAGE REPEAT', 'SAG MP-RAGE'
                                                       'MT1; GradWarp; N3m <- IR-FSPGR', ]
                                )

                        descriptions = set(descriptions)
                        is_accept = list(accepted_descriptions.intersection(descriptions))

                        if not is_accept:
                            # Ask the user to choose one
                            print('Please choose one:')
                            mri_image_id = input('Enter the image ID for subject %s: ' % subject)
                            print('You entered: %s' % mri_image_id)

                            selected = age_df.loc[age_df['Image ID'] == int(mri_image_id)]
                            id = selected['Image ID'].values[0]
                            description = selected['Description'].values[0]
                            new_row = pd.DataFrame(
                                    [{'SUBJECT_ID': subject, 'MRI_AGE': age, 'MRI_ID': id, 'DESCRIPTION': description}],
                                    columns = ['SUBJECT_ID', 'MRI_AGE', 'MRI_ID', 'DESCRIPTION']
                                    )
                            filtered_mri_df = pd.concat([filtered_mri_df, new_row], ignore_index = True)

                        else:
                            selected = age_df.loc[age_df['Description'] == is_accept[0]]
                            print(selected)
                            id = selected['Image ID'].values[0]
                            description = selected['Description'].values[0]
                            new_row = pd.DataFrame(
                                    [{'SUBJECT_ID': subject, 'MRI_AGE': age, 'MRI_ID': id, 'DESCRIPTION': description}],
                                    columns = ['SUBJECT_ID', 'MRI_AGE', 'MRI_ID', 'DESCRIPTION']
                                    )
                            filtered_mri_df = pd.concat([filtered_mri_df, new_row], ignore_index = True)

                        print('\n')

        filtered_mri_df.to_csv(out, index = False)

    def is_subject_amyloid_PET_positive(self, subject_id = False, rid = False, date = 'baseline'):
        """
        Checks if a subject is amyloid PET positive.

        :param subject_id: The ID of the subject.
        :param rid: The RID of the subject.
        :param date: The date of the PET scan.
        :return: The tracer, SUVR, amyloid status, and centiloids.
        """

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
        """
        Gets the CSF biomarkers for a subject.

        :param subject_id: The ID of the subject.
        :param rid: The RID of the subject.
        :param date: The date of the CSF biomarker measurement (The measurement closest to this date will be returned.
        :return: The AB42, tau, ptau, asyn, and exam date.
        """
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

            ordered = df_subj.sort_values(by = ['EXAMDATE'], ascending = True)

            if date == 'baseline':
                ab42 = ordered['ABETA42'].values[0]
                tau = ordered['TAU'].values[0]
                ptau = ordered['PTAU'].values[0]

            elif date == 'last':
                ab42 = ordered['ABETA42'].values[-1]
                tau = ordered['TAU'].values[-1]
                ptau = ordered['PTAU'].values[-1]
            else:
                date = datetime.strptime(date, '%Y-%m-%d')
                diff_time_csf = 1000

                for csf_index, csf_row in ordered.iterrows():

                    date_csf = ordered["EXAMDATE"]
                    ab42_ = ordered["ABETA42"]
                    tau_ = ordered["TAU"]
                    ptau_ = ordered["PTAU"]

                    diff_time_ = (
                            relativedelta(date_csf, date).years
                            + relativedelta(date_csf, date).months / 12
                            + relativedelta(date_csf, date).days / 365
                    )

                    if np.abs(diff_time_) < diff_time_csf:
                        diff_time_csf = np.abs(diff_time_)
                        ab42 = ab42_
                        tau = tau_
                        ptau = ptau_
                        exam_date = date_csf

        df_subj = self.asyn_csf_df.loc[self.asyn_csf_df['RID'] == rid]

        if not df_subj.empty:

            ordered = df_subj.sort_values(by = ['EXAMDATE'], ascending = True)

            if date == 'baseline':
                asyn_in = ordered['Result'].values[0]
                exam_date = ordered['EXAMDATE'].values[0]
            elif date == 'last':
                asyn_in = ordered['Result'].values[-1]
                exam_date = ordered['EXAMDATE'].values[-1]
            else:
                date = datetime.strptime(date, '%Y-%m-%d')
                diff_time_csf = 1000

                for csf_index, csf_row in ordered.iterrows():

                    date_csf = ordered["EXAMDATE"]
                    asyn_ = ordered["Result"]

                    diff_time_ = (
                            relativedelta(date_csf, date).years
                            + relativedelta(date_csf, date).months / 12
                            + relativedelta(date_csf, date).days / 365
                    )

                    if np.abs(diff_time_) < diff_time_mmse:
                        diff_time_mmse = np.abs(diff_time_)
                        asyn_in = asyn_


            if asyn_in == 'Not_Detected':

                asyn = 0

            elif asyn_in == 'Detected-1':

                asyn = 1

            elif asyn_in == 'Detected-2':

                asyn = 2

            else:
                asyn = None

        return ab42, tau, ptau, asyn, exam_date

    def get_genetics_data(self, subject_id = False, rid = False):
        """
        Gets the genetic data for a subject.

        :param subject_id: The ID of the subject.
        :param rid: The RID of the subject.
        :return: The APOE, KCNMB2, TMEM106B, GRN, rs73069071, and ABCC9.
        """

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

        if not filtered_mmse.empty:

            ordered = filtered_mmse.sort_values(by = ['USERDATE'], ascending = True)

            if date == 'baseline':
                mmse = ordered['MMSCORE'].values[0]
            elif date == 'last':
                mmse = ordered['MMSCORE'].values[-1]
            else:
                date = datetime.strptime(date, '%Y-%m-%d')
                diff_time_mmse = 1000

                for mmse_index, mmse_row in ordered.iterrows():

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

        filtered_composite = self.composites_df[self.composites_df.RID == rid]
        filtered_composite['EXAMDATE'] = pd.to_datetime(filtered_composite.EXAMDATE)

        if not filtered_composite.empty:

            ordered = filtered_composite.sort_values(by = ['EXAMDATE'], ascending = True)

            if date == 'baseline':
                adni_mem = ordered['ADNI_MEM'].dropna().values[0]
                adni_exec = ordered['ADNI_EF'].dropna().values[0]
                adni_lan = ordered['ADNI_LAN'].dropna().values[0]
                adni_vs = ordered['ADNI_VS'].dropna().values[0]
            elif date == 'last':
                adni_mem = ordered['ADNI_MEM'].dropna().values[-1]
                adni_exec = ordered['ADNI_EF'].dropna().values[-1]
                adni_lan = ordered['ADNI_LAN'].dropna().values[-1]
                adni_vs = ordered['ADNI_VS'].dropna().values[-1]
            else:
                #date = datetime.strptime(date, '%Y-%m-%d')
                diff_time_composite = 1000

                for composite_index, composite_row in ordered.iterrows():

                    date_composite = composite_row["EXAMDATE"]
                    adni_mem_ = composite_row["ADNI_MEM"]
                    adni_exec_ = composite_row["ADNI_EF"]
                    adni_lan_ = composite_row["ADNI_LAN"]
                    adni_vs_ = composite_row["ADNI_VS"]

                    diff_time_ = (
                            relativedelta(date_composite, date).years
                            + relativedelta(date_composite, date).months / 12
                            + relativedelta(date_composite, date).days / 365
                    )

                    if np.abs(diff_time_) < diff_time_composite:
                        diff_time_composite = np.abs(diff_time_)
                        adni_mem = adni_mem_
                        adni_exec = adni_exec_
                        adni_lan = adni_lan_
                        adni_vs = adni_vs_

        return mmse, adni_mem, adni_exec, adni_lan, adni_vs

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