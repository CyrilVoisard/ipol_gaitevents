# Objectif : construire une classe qui permet de sélectionner des fichiers

import os
import numpy as np
import pickle

home = os.getcwd()


class Patient(object):

    def __init__(self, id, home):
        self.id = id
        self.name = id[4:10]
        self.cohorte = "PVE"
        self.path = os.path.join(home, id)

        # Liste des mois
        self.liste_mois = list()

        # Liste des dates
        self.dates = []

        # Liste des fichiers
        self.gr_exp = []
        self.xs_exp = []

        if os.path.exists(os.path.join(self.path, "M0")):
            # Fichiers de signaux
            self.liste_mois = self.liste_mois + ["M0"]
            self.path_M0_GR = os.path.join(self.path, "M0", "GR")
            for file in os.listdir(self.path_M0_GR):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M0", "GR", file[:26]))
            self.path_M0_XS = os.path.join(self.path, "M0", "XSens")
            for file in os.listdir(self.path_M0_XS):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M0", "XSens", file[:26]))
                    else:
                        if file[20] == 'X':
                            self.xs_exp.append(os.path.join("M0", "XSens", file[:26]))
            # Liste de dates
            for file in os.listdir(self.path_M0_XS):
                if file[0:6] == self.name:
                    date = file[11:19]
                    if date not in self.dates:
                        self.dates.append(date)
        if os.path.exists(os.path.join(self.path, "M3")):
            self.liste_mois = self.liste_mois + ["M3"]
            self.path_M3_GR = os.path.join(self.path, "M3", "GR")
            for file in os.listdir(self.path_M3_GR):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M3", "GR", file[:26]))
            self.path_M3_XS = os.path.join(self.path, "M3", "XSens")
            for file in os.listdir(self.path_M3_XS):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M3", "XSens", file[:26]))
                    else:
                        if file[20] == 'X':
                            self.xs_exp.append(os.path.join("M3", "XSens", file[:26]))
            # Liste de dates
            for file in os.listdir(self.path_M3_XS):
                if file[0:6] == self.name:
                    date = file[11:19]
                    if date not in self.dates:
                        self.dates.append(date)
        if os.path.exists(os.path.join(self.path, "M6")):
            self.liste_mois = self.liste_mois + ["M6"]
            self.path_M6_GR = os.path.join(self.path, "M6", "GR")
            for file in os.listdir(self.path_M6_GR):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M6", "GR", file[:26]))
            self.path_M6_XS = os.path.join(self.path, "M6", "XSens")
            for file in os.listdir(self.path_M6_XS):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M6", "XSens", file[:26]))
                    else:
                        if file[20] == 'X':
                            self.xs_exp.append(os.path.join("M6", "XSens", file[:26]))
            # Liste de dates
            for file in os.listdir(self.path_M6_XS):
                if file[0:6] == self.name:
                    date = file[11:19]
                    if date not in self.dates:
                        self.dates.append(date)
        if os.path.exists(os.path.join(self.path, "M9")):
            self.liste_mois = self.liste_mois + ["M9"]
            self.path_M9_GR = os.path.join(self.path, "M9", "GR")
            for file in os.listdir(self.path_M9_GR):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M9", "GR", file[:26]))
            self.path_M9_XS = os.path.join(self.path, "M9", "XSens")
            for file in os.listdir(self.path_M9_XS):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(os.path.join("M9", "XSens", file[:26]))
                    else:
                        if file[20] == 'X':
                            self.xs_exp.append(os.path.join("M9", "XSens", file[:26]))
            # Liste de dates
            for file in os.listdir(self.path_M9_XS):
                if file[0:6] == self.name:
                    date = file[11:19]
                    if date not in self.dates:
                        self.dates.append(date)

        # Détermination des dates de prise
        self.dates.sort()
        #print("Le patient a fait ces enregistrements : ", self.liste_mois)
        #print("Le patient a fait des enregistrements à ces dates : ", self.dates)

        #print(self.gr_exp)
        self.gr_exp = list(set(self.gr_exp))
        self.gr_exp.sort()
        self.gr_exp = np.array(self.gr_exp)
        self.xs_exp = list(set(self.xs_exp))
        self.xs_exp.sort()
        self.xs_exp = np.array(self.xs_exp)

        self.gr_ref = self.gr_exp
        # self.gr_ref = self.gr_exp[self.gr_exp[:, 0] == id + "_" + self.dates[0]]

        # Existence d'un dictionnaire de pas
        if (os.path.exists(os.path.join(self.path, 'verif_gr'))):
            liste_racine = os.listdir(os.path.join(self.path, 'verif_gr'))
            liste_dict = [l for l in liste_racine if l.endswith("dictionnaire.pkl")]
            if len(liste_dict) == 1 :
                self.is_dictionnaire = True
                # Ouverture du dictionnaire
                dic_path = os.path.join(self.path, 'verif_gr', liste_dict[0])
                with open(dic_path, 'rb') as f:
                    dict_open = pickle.load(f)
                self.dictionnaire = dict_open
            else:
                self.is_dictionnaire = False
        else:
            self.is_dictionnaire = False

        self.meta = dict()
