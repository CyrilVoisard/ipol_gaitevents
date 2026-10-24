import os
import json
from load_package import xsens, gait_rite


class Experience_GR(object):

    # ---------------------------------- Initialisation ---------------------------------------------

    def __init__(self, root, id_exp, freq=100):

        # Racine
        self.root = root
        # Nom
        self.name = id_exp[0:6]
        # Pathologie
        self.patho = id_exp[7:10]
        # Date
        self.date = id_exp[11:19]
        # Type d'expérience
        self.type = id_exp[20:26]
        # Id_exp
        self.id_exp = id_exp
        # Fréquence d'acquisition
        self.freq = freq

        # Liste des fichiers de l'expérience
        root_list = os.listdir(root)
        # print(root_list, id_exp)
        self.file_list = [e for e in root_list if e.startswith(id_exp)]

        # Existe-t-il une segmentation des pas
        self.steps_list = [s for s in self.file_list if s.endswith('_steps.txt')]
        self.steps_list_auto = [s for s in self.steps_list if s.endswith('auto_steps.txt')]
        self.steps_list_manu = [s for s in self.steps_list if s.endswith('manu_steps.txt')]
        self.steps_list_alnr = [s for s in self.steps_list if s.endswith('alnr_steps.txt')]
        self.steps_list_semi = [s for s in self.steps_list if s.endswith('semi_steps.txt')]

        self.is_steps = []
        self.path_steps = dict()

        if len(self.steps_list_auto) + len(self.steps_list_manu) + len(self.steps_list_alnr) + len(
                self.steps_list_semi) == 0:
            self.is_steps.append('aucun')
        else:
            if len(self.steps_list_auto) > 0:
                self.is_steps.append('auto')
                self.path_steps["auto"] = os.path.join(self.root, self.steps_list_auto[0])
            if len(self.steps_list_manu) > 0:
                self.is_steps.append('manu')
                self.path_steps["manu"] = os.path.join(self.root, self.steps_list_manu[0])
            if len(self.steps_list_semi) > 0:
                self.is_steps.append('semi')
                self.path_steps["semi"] = os.path.join(self.root, self.steps_list_semi[0])
            if len(self.steps_list_alnr) > 0:
                self.is_steps.append('alnr')
                self.path_steps["alnr"] = os.path.join(self.root, self.steps_list_alnr[0])

        # Existe-t-il une segmentation des phases de la marche
        self.seg_list = [s for s in self.file_list if s.endswith('_seg.txt')]
        self.seg_list_auto = [s for s in self.seg_list if s.endswith('auto_seg.txt')]
        self.seg_list_manu = [s for s in self.seg_list if s.endswith('manu_seg.txt')]

        self.is_seg = []

        if len(self.seg_list_auto) + len(self.seg_list_manu) == 0:
            self.is_seg.append('aucun')
        else:
            if len(self.seg_list_auto) > 0:
                self.is_seg.append('auto')
            if len(self.seg_list_manu) > 0:
                self.is_seg.append('manu')

    # ---------------------------------- Importations ---------------------------------------------

    def get_XSens_brut(self, capteur, disp=False):
        file_capteur = [d for d in self.file_list if d.endswith(capteur)][0]
        path_capteur = os.path.join(self.root, file_capteur)
        data = xsens.load_XSens_brut(path_capteur)
        return data

    def get_data_tronc_brut(self, disp=False):

        with open("sensors_dictionnaires.json") as json_dict:
            sensors_dict = json.load(json_dict)

        i = 0
        j = 0

        while (i < len(sensors_dict)) & (j != 10):
            try:
                data_tronc = self.get_XSens_brut(sensors_dict[i]['tronc'], disp)
                j = 10
            except:
                i = i + 1
        if j == 0:
            print(self.id_exp)

        return data_tronc

    def get_XSens(self, capteur, t0_min=float('inf'), correction_acc=True, start=0, end=200, filtre=True, N=8, fc=14,
                  disp=False):
        file_capteur = [d for d in self.file_list if d.endswith(capteur)][0]
        path_capteur = os.path.join(self.root, file_capteur)
        # print(capteur, path_capteur)
        data = xsens.import_XSens(path_capteur, t0_min, correction_acc, start, end, filtre, N, fc, freq=self.freq,
                                  disp=disp)
        return data

    def get_data_tronc(self, correction_acc=True, disp=False):

        with open("sensors_dictionnaires.json") as json_dict:
            sensors_dict = json.load(json_dict)

        i = 0
        j = 0
        while (i < len(sensors_dict)) & (j != 10):
            try:
                data_tronc = self.get_XSens(sensors_dict[i]['tronc'], correction_acc=correction_acc, disp=disp)
                j = 10
            except:
                i = i + 1
        if j == 0:
            print(self.id_exp)

        return data_tronc

    def import_experience(self, correction_acc=True, correction_gr=True, disp=False):

        with open("sensors_dictionnaires.json") as json_dict:
            sensors_dict = json.load(json_dict)

        i = 0
        j = 0

        while (i < len(sensors_dict)) & (j != 10):
            try:
                data_tete_brut = self.get_XSens_brut(sensors_dict[i]['tete'])
                data_tronc_brut = self.get_XSens_brut(sensors_dict[i]['tronc'])
                data_pied_droit_brut = self.get_XSens_brut(sensors_dict[i]['pied_droit'])
                data_pied_gauche_brut = self.get_XSens_brut(sensors_dict[i]['pied_gauche'])
                t0_min = xsens.find_t0_min(
                    [data_tete_brut, data_tronc_brut, data_pied_droit_brut, data_pied_gauche_brut])
                data_tete = self.get_XSens(sensors_dict[i]['tete'], t0_min, correction_acc=correction_acc, disp=disp)
                data_tronc = self.get_XSens(sensors_dict[i]['tronc'], t0_min, correction_acc=correction_acc, disp=disp)
                data_pied_droit = self.get_XSens(sensors_dict[i]['pied_droit'], t0_min, correction_acc=correction_acc,
                                                 disp=disp)
                data_pied_gauche = self.get_XSens(sensors_dict[i]['pied_gauche'], t0_min, correction_acc=correction_acc,
                                                  disp=disp)
                j = 10
            except:
                i = i + 1
        if j == 0:
            print(self.id_exp)

        # print(self.file_list)
        if not correction_gr:
            file_gait = [d for d in self.file_list if d.endswith("CoP.txt")][0]
            path_gait = os.path.join(self.root, file_gait)
            exo = gait_rite.process_gaitrite(path_gait)
        else:
            try:
                file_gait = [d for d in self.file_list if d.endswith("CoP_relu.txt")][0]
                path_gait = os.path.join(self.root, file_gait)
                exo = gait_rite.process_gaitrite_corrected(path_gait)
            except:
                print("Pas possible de corriger le GR !")

        return data_tete, data_tronc, data_pied_droit, data_pied_gauche, exo

    def import_experience_pied(self, correction_acc=True, correction_gr=True, disp=False):
        #print("Correction : ", correction_acc)

        with open("sensors_dictionnaires.json") as json_dict:
            sensors_dict = json.load(json_dict)

        i = 0
        j = 0

        while (i < len(sensors_dict)) & (j != 10):
            try:
                # print(sensors_dict[i]['pied_droit'])
                data_pied_droit_brut = self.get_XSens_brut(sensors_dict[i]['pied_droit'])
                data_pied_gauche_brut = self.get_XSens_brut(sensors_dict[i]['pied_gauche'])
                t0_min = xsens.find_t0_min([data_pied_droit_brut, data_pied_gauche_brut])
                data_pied_droit = self.get_XSens(sensors_dict[i]['pied_droit'], t0_min, correction_acc=correction_acc,
                                                 disp=disp)
                data_pied_gauche = self.get_XSens(sensors_dict[i]['pied_gauche'], t0_min, correction_acc=correction_acc,
                                                  disp=disp)
                # print("ok")
                j = 10
            except:
                i = i + 1
        if j == 0:
            print(self.id_exp)

        # print(self.file_list)

        if not correction_gr:
            file_gait = [d for d in self.file_list if d.endswith("CoP.txt")][0]
            path_gait = os.path.join(self.root, file_gait)
            exo = gait_rite.process_gaitrite(path_gait)
        else:
            # try:
            file_gait = [d for d in self.file_list if d.endswith("CoP_relu.txt")][0]
            # print(file_gait)
            path_gait = os.path.join(self.root, file_gait)
            # print(os.path.exists(path_gait), path_gait)
            exo = gait_rite.process_gaitrite_corrected(path_gait)
            # except:
            #   print("Pas possible de corriger le GR !")

        return data_pied_droit, data_pied_gauche, exo

    def import_exo(self, correction_gr=True, disp=False):
        if not correction_gr:
            file_gait = [d for d in self.file_list if d.endswith("CoP.txt")][0]
            path_gait = os.path.join(self.root, file_gait)
            exo = gait_rite.process_gaitrite(path_gait)
        else:
            try:
                file_gait = [d for d in self.file_list if d.endswith("CoP_relu.txt")][0]
                path_gait = os.path.join(self.root, file_gait)
                exo = gait_rite.process_gaitrite_corrected(path_gait)
            except:
                print("Pas possible de corriger le GR !")
        return exo

    def import_tapis(self, disp=False):
        file_gait = [d for d in self.file_list if d.endswith("CoP.txt")][0]
        path_gait = os.path.join(self.root, file_gait)
        tapis = gait_rite.importfile_GaitRite(path_gait)
        return tapis

    def get_steps(self, steps_pref='semi', disp=False):  # Segmentation des pas : automatique meilleure que manuelle

        steps_lim = 0

        if len(self.is_steps) == 0:
            if disp:
                print('Problème de fichier, valeur retournée : 0.')

        else:
            try:
                steps_path = self.path_steps[steps_pref]
                steps_lim = seg_steps.import_gait_events(steps_path, disp)
            except:
                if disp:
                    print('Segmentation demandée non disponible, première segmentation par défaut dans la liste : ',
                          self.is_steps)
                steps_path = list(self.path_steps.values())[0]
                steps_lim = seg_steps.import_gait_events(steps_path, disp)

        return steps_lim

    def get_seg(self, seg_pref='manu', disp=False):
        # Segmentation des pas : manuelle meilleure que automatique si pas de segmentation des pas
        # (pour le début et la fin)

        seg_lim = 0

        if self.is_seg in ['probleme', 'aucun']:
            if disp:
                print('Problème de fichier, valeur retournée : 0.')

        else:
            if seg_pref == 'manu':
                try:
                    seg_path = os.path.join(self.root, self.seg_list_manu[0])
                    seg_lim = seg_steps.import_seg(seg_path, disp)
                except:
                    if disp:
                        print('Pas de segmentation manuelle, segmentation automatique par défaut.')
                    seg_path = os.path.join(self.root, self.seg_list_auto[0])
                    seg_lim = seg_steps.import_seg(seg_path, disp)
            else:
                if seg_pref == 'auto':
                    try:
                        seg_path = os.path.join(self.root, self.seg_list_auto[0])
                        seg_lim = seg_steps.import_seg(seg_path, disp)
                    except:
                        if disp:
                            print('Pas de segmentation automatique, segmentation manuelle par défaut.')
                        seg_path = os.path.join(self.root, self.seg_list_manu[0])
                        seg_lim = seg_steps.import_seg(seg_path, disp)

        return seg_lim
