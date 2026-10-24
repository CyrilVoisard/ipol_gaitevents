# Objectif : construire une classe qui permet de sélectionner des fichiers

import os
import numpy as np
import pickle


class Patient(object):

    def __init__(self, id, home):
        self.id = id
        self.name = id[0:6]
        self.cohorte = id[7:10]
        self.path = os.path.join(home, id)

        # Détermination des dates de prise
        self.dates = []
        for file in os.listdir(self.path):
            if file[0:6] == self.name:
                date = file[11:19]
                if date not in self.dates:
                    self.dates.append(date)
        self.dates.sort()

        # Liste des fichiers
        self.gr_exp = []
        self.xs_exp = []

        for file in os.listdir(self.path):
            if file.endswith(".txt"):
                if len(file) > 21:
                    if file[20] == 'G':
                        self.gr_exp.append(file[0:26])
                    else:
                        if (file[20] == 'X') | (file[20] == 'P'):
                            self.xs_exp.append(file[0:26])
                        else:
                            self.xs_exp.append(file[0:-11])

        self.gr_exp = np.array(self.gr_exp)
        self.gr_exp = list(set(self.gr_exp))
        self.gr_exp.sort()
        self.gr_exp = np.array(self.gr_exp)
        self.xs_exp = list(set(self.xs_exp))
        self.xs_exp.sort()
        self.xs_exp = np.array(self.xs_exp)

        # if len(self.gr_exp) > 0:
        # self.gr_ref = self.gr_exp[self.gr_exp[:, 0] == id + "_" + self.dates[0]]

        # Existence d'un dictionnaire de pas
        if os.path.exists(os.path.join(self.path, 'verif_gr')):
            liste_racine = os.listdir(os.path.join(self.path, 'verif_gr'))
            liste_dict = [l for l in liste_racine if l.endswith("dictionnaire.pkl")]
            if len(liste_dict) == 1:
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
