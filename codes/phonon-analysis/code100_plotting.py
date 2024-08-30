# KG 7/11/2021

import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from code2_ev_and_freq_from_phonopy import phononBC
from code2_ev_and_freq_for_gamma import phononBC as phononBC_for_gamma
from code4_modal_calc_for_non_gamma import En_hist_n_center

class contributing_phonons_n_PBC_NEB:
    def __init__(self,val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15, val16, val17, val18, val19):
        self.w = val1 # Phonon frequencies
        self.v = val2 # Phonon eigenvectors
        self.En = val3 # modal contributions
        self.atomsp = val4 # atom species (string)
        self.atomNs = val5 # atom species numbers (int)
        self.compound_name = val6 # compound_name (for plot saving)
        self.bn = val7 # modal bottleneckk area
        self.PBCTotal = val8
        self.PBCO = val9
        self.PBCA = val10
        self.PBCB = val11
        self.counter_TotalDOS = val12
        self.counter_partialDOS_O = val13
        self.counter_partialDOS_A = val14
        self.counter_partialDOS_B = val15
        self.w_TotalDOS = val16
        self.w_partialDOS_O = val17
        self.w_partialDOS_A = val18
        self.w_partialDOS_B = val19

    def plot_En_PBC(self, AB_first,filename1):
        freq = self.w.copy()
        self.En /= np.sum(self.En) # Normalizing En
        En = self.En.copy()
        TotalPBC = self.PBCTotal.copy()
        partialPBC = self.PBCO.copy()
        partialPBC_A = self.PBCA.copy()
        partialPBC_B = self.PBCB.copy()
        counter_TotalDOS = self.counter_TotalDOS.copy()
        counter_partialDOS_O = self.counter_partialDOS_O.copy()
        counter_partialDOS_A = self.counter_partialDOS_A.copy()
        counter_partialDOS_B = self.counter_partialDOS_B.copy()
        w_TotalDOS = self.w_TotalDOS.copy()
        w_partialDOS_O = self.w_partialDOS_O.copy()
        w_partialDOS_A = self.w_partialDOS_A.copy()
        w_partialDOS_B = self.w_partialDOS_B.copy()

        modalEn_instance = En_hist_n_center(freq, En, 50)
        [w_En, counter_En, w_En_center] = modalEn_instance.hist_for_En_n_center()
        #phononBC_instance = phononBC(50, freq, self.v, self.atomNs, self.atomsp)
        # [w_TotalDOS, counter_TotalDOS, TotalPBC] = phononBC_instance.TotalDOS() # Calculate Total DOS histogram
        # [w_partialDOS, counter_partialDOS_O, partialPBC] = phononBC_instance.partialDOS('O', AB_first) # Calculate partial DOS histogram
        # [w_partialDOS_A, counter_partialDOS_A, partialPBC_A] = phononBC_instance.partialDOS('A', AB_first)
        # [w_partialDOS_B, counter_partialDOS_B, partialPBC_B] = phononBC_instance.partialDOS('B', AB_first)
        counter_partialDOS_O /= np.sum(counter_TotalDOS) # Normalizing partial DOS first
        counter_partialDOS_A /= np.sum(counter_TotalDOS)
        counter_TotalDOS /= np.sum(counter_TotalDOS) # Then, normalizing total DOS
        
        # plot series 1
        print(np.sum(counter_partialDOS_A))
        plt.figure(self.compound_name)
        plt.rcParams.update({'font.size': 13})
        plt.fill_between(w_En, 0, self.smooth(counter_En,3), color="black", alpha=0.5, label='Phonon contributions')
        plt.fill_between(w_TotalDOS, 0, self.smooth(counter_TotalDOS,3), color="blue", alpha=0.5, label='Total phonon DOS')
        plt.fill_between(w_partialDOS_O, 0, self.smooth(counter_partialDOS_O,3), color="red", alpha=0.5, label='Anion Partial pDOS')
        plt.fill_between(w_partialDOS_A, 0, self.smooth(counter_partialDOS_A,3), color="green", alpha=0.5, label='Cation Partial pDOS')
        plt.plot(freq,En, linestyle='none', marker='o', markerfacecolor='none', markersize=7, markeredgecolor='k', label='Phonon contributions')

        filename = filename1+"phonon_total_DOS_" + self.compound_name + ".txt"
        f_total_DOS = open(filename, "w")
        for i in np.arange(0,w_TotalDOS.shape[0],1):
            f_total_DOS.write("%f %f\n" %(w_TotalDOS[i], self.smooth(counter_TotalDOS,3)[i]))
        f_total_DOS.close()

        filename = filename1+"phonon_partial_DOS_" + self.compound_name + ".txt"
        f_partial_DOS = open(filename, "w")
        for i in np.arange(0,w_partialDOS_O.shape[0],1):
            f_partial_DOS.write("%f %f\n" %(w_partialDOS_O[i], self.smooth(counter_partialDOS_O,3)[i]))
        f_partial_DOS.close()

        minx = np.floor(np.min(freq))
        maxx = np.ceil(np.max(freq))
        miny = 0
        maxy = np.ceil(np.max(En)*20.)/20.
        print(np.max(En))

        plt.plot([w_En_center,w_En_center],[miny,0.6*maxy],'k--')
        plt.plot([TotalPBC,TotalPBC],[miny,0.6*maxy],'b--')
        plt.plot([partialPBC,partialPBC],[miny,0.6*maxy],'r--')
        plt.plot([partialPBC_A,partialPBC_A],[miny,0.6*maxy],'g--')
        plt.text(w_En_center, 0.63*maxy, 'En_center', fontsize=10, rotation = 'vertical', color = 'k')
        plt.text(TotalPBC, 0.63*maxy, 'tPBC', fontsize=10, rotation = 'vertical', color = 'b')
        plt.text(partialPBC, 0.63*maxy, 'Anion PBC', fontsize=10, rotation = 'vertical', color = 'r')
        plt.text(partialPBC_A, 0.63*maxy, 'Cation PBC', fontsize=10, rotation = 'vertical', color = 'g')
        plt.xlabel('Frequency (cm$^{-1}$)')
        plt.ylabel('Normlized excitation energy')
        plt.xlim(minx, maxx)
        plt.ylim(miny, maxy)

        plt.legend(fontsize=11)
        # filename1 = 'phonon_' + self.compound_name + '.png'
        plt.savefig(filename1+'phonon.png', bbox_inches = 'tight')
        #plt.show()
        plt.close()

        
                ##### Saving phonon contribution (En) data
        # En
        column_stacked_variables = np.column_stack((freq, En))
        filename2 =filename1+ 'phonon_contributions_to_hop.csv'
        np.savetxt(filename2, column_stacked_variables, delimiter=",")

        # En DOS
        column_stacked_variables = np.column_stack((w_En, counter_En))
        filename2 = filename1+'phonon_contributions_DOS_to_hop.csv'
        np.savetxt(filename2, column_stacked_variables, delimiter=",")
        column_stacked_variables = [w_En_center,TotalPBC,partialPBC,partialPBC_A]
        filename2 = filename1+'phonon_contributions_DOS_to_hop_center.csv'
        np.savetxt(filename2, column_stacked_variables, delimiter=",")
        # stuff related to plotting bottleneck area change
#         self.bn /= np.sum(self.bn) # Normalizing bn
#         bn = self.bn.copy()
#         modalbn_instance = En_hist_n_center(freq, bn, 50)
#         [w_bn, counter_bn, w_bn_center] = modalbn_instance.hist_for_En_n_center()

#         # plot series 2
#         plt.figure(self.compound_name)
#         plt.rcParams.update({'font.size': 13})
#         plt.fill_between(w_En, 0, self.smooth(counter_En,3), color="black", alpha=0.5, label='Phonon contributions')
#         plt.plot(freq,En, linestyle='none', marker='o', markerfacecolor='none', markersize=7, markeredgecolor='k', label='Phonon contributions')
#         plt.fill_between(w_bn, 0, self.smooth(counter_bn,3), color="red", alpha=0.5, label='modal change in bottleneck')
#         plt.plot(freq,bn, linestyle='none', marker='s', markerfacecolor='none', markersize=7, markeredgecolor='r', label='modal change in bottleneck')

#         minx = np.floor(np.min(freq))
#         maxx = np.ceil(np.max(freq))
#         miny = np.floor(np.min(bn)*20.)/20.
#         maxy = np.ceil(np.max(En)*20.)/20.

#         plt.plot([w_En_center,w_En_center],[0,0.6*maxy],'k--')
#         plt.plot([w_bn_center,w_bn_center],[0,0.6*maxy],'r--')
#         #plt.text(w_En_center, 0.63*maxy, 'En_center', fontsize=10, rotation = 'vertical', color = 'k')
#         #plt.text(w_bn_center, 0.63*maxy, 'bn_change_center', fontsize=10, rotation = 'vertical', color = 'k')

#         plt.xlabel('Frequency (cm$^{-1}$)')
#         plt.ylabel('Normlized excitation energy and\nbottleneck area change')
#         plt.xlim(minx, maxx)
#         plt.ylim(miny, maxy)

#         # Make room for the labels
#         plt.gcf().subplots_adjust(bottom=0.15)
#         plt.gcf().subplots_adjust(left=0.2)

#         plt.legend(fontsize=11)
#         filename = 'area_' + self.compound_name + '.png'
#         #plt.savefig(filename, bbox_inches = 'tight')
#         #plt.show()
#         plt.close()

        # plot series 3
#         plt.figure(self.compound_name)
#         plt.rcParams.update({'font.size': 13})
#         plt.fill_between(w_bn, 0, self.smooth(counter_bn,3), color="red", alpha=0.5, label='modal change in bottleneck')
#         plt.plot(freq,bn, linestyle='none', marker='s', markerfacecolor='none', markersize=7, markeredgecolor='r', label='modal change in bottleneck')

#         minx = np.floor(np.min(freq))
#         maxx = np.ceil(np.max(freq))
#         miny = np.floor(np.min(bn)*100.)/100.
#         maxy = np.ceil(np.max(bn)*12.)/12.

#         #plt.plot([w_bn_center,w_bn_center],[0,0.6*maxy],'r--')
#         #plt.text(w_En_center, 0.63*maxy, 'En_center', fontsize=10, rotation = 'vertical', color = 'k')
#         #plt.text(w_bn_center, 0.63*maxy, 'bn_change_center', fontsize=10, rotation = 'vertical', color = 'k')

#         plt.xlabel('Frequency (cm$^{-1}$)')
#         plt.ylabel('Normlized bottleneck area change\nby each phonon')
#         plt.xlim(minx, maxx)
#         plt.ylim(miny, maxy)

#         # Make room for the labels
#         plt.gcf().subplots_adjust(bottom=0.15)
#         plt.gcf().subplots_adjust(left=0.2)

#         plt.legend(fontsize=11)
#         filename = 'pure_area_' + self.compound_name + '.png'
#         #plt.savefig(filename, bbox_inches = 'tight')
#         #plt.show()
#         plt.close()

        # plot series 4
        plt.figure(self.compound_name)
        plt.rcParams.update({'font.size': 13})
        #plt.fill_between(w_En, 0, self.smooth(counter_En,3), color="black", alpha=0.5, label='Phonon contributions')
        #plt.fill_between(w_TotalDOS, 0, self.smooth(counter_TotalDOS,3), color="blue", alpha=0.5, label='Total phonon DOS')
        #plt.fill_between(w_partialDOS, 0, self.smooth(counter_partialDOS,3), color="red", alpha=0.5, label='Oxygen pDOS')
        #plt.plot(freq,En, linestyle='none', marker='o', markerfacecolor='none', markersize=7, markeredgecolor='k', label='Phonon contributions')

        minx = np.floor(np.min(freq))
        maxx = 50
        miny = 0
        maxy = np.ceil(np.max(counter_TotalDOS)*20.)/20.
        print(np.max(counter_TotalDOS))

        #plt.plot([w_En_center,w_En_center],[miny,0.6*maxy],'k--')
        plt.plot([TotalPBC,TotalPBC],[miny,0.6*maxy],'b--')
        plt.plot([partialPBC,partialPBC],[miny,0.6*maxy],'r--')

        #plt.text(w_En_center, 0.63*maxy, 'En_center', fontsize=10, rotation = 'vertical', color = 'k')
        # plt.text(TotalPBC, 0.63*maxy, 'tPBC', fontsize=10, rotation = 'vertical', color = 'b')
        # plt.text(partialPBC, 0.63*maxy, 'pPBC', fontsize=10, rotation = 'vertical', color = 'r')
        
        plt.xlabel('Frequency (cm$^{-1}$)')
        plt.ylabel('Phonon DOS')
        plt.xlim(minx, maxx)
        plt.ylim(miny, maxy)

        #plt.legend(fontsize=11)
        #filename = 'PDOS_' + self.compound_name + '.png'
        #plt.savefig(filename, bbox_inches = 'tight')
        ##plt.show()
        plt.close()

        '''
        ##### Saving PDOS data in different files
        # Total PDOS
        column_stacked_variables = np.column_stack((w_TotalDOS, counter_TotalDOS))
        filename = 'PDOS_' + self.compound_name + '_total.csv'
        np.savetxt(filename, column_stacked_variables, delimiter=",")
        # Partial PDOS O 
        column_stacked_variables = np.column_stack((w_partialDOS, counter_partialDOS))
        filename = 'PDOS_' + self.compound_name + '_partial_O.csv'
        np.savetxt(filename, column_stacked_variables, delimiter=",")
        # Partial PDOS A
        column_stacked_variables = np.column_stack((w_partialDOS_A, counter_partialDOS_A))
        filename = 'PDOS_' + self.compound_name + '_partial_A.csv'
        np.savetxt(filename, column_stacked_variables, delimiter=",")
        # Partial PDOS B 
        column_stacked_variables = np.column_stack((w_partialDOS_B, counter_partialDOS_B))
        filename = 'PDOS_' + self.compound_name + '_partial_B.csv'
        np.savetxt(filename, column_stacked_variables, delimiter=",")
        
        ##### Saving phonon contribution (En) data
        # En
        column_stacked_variables = np.column_stack((freq, En))
        filename = 'phonon_contributions_to_hop_' + self.compound_name + '.csv'
        np.savetxt(filename, column_stacked_variables, delimiter=",")

        # En DOS
        column_stacked_variables = np.column_stack((w_En, counter_En))
        filename = 'phonon_contributions_DOS_to_hop_' + self.compound_name + '.csv'
        np.savetxt(filename, column_stacked_variables, delimiter=",")
        '''
        
    def smooth(self, y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

class contributing_phonons_n_PBC_non_NEB:
    def __init__(self,val1, val2, val4, val5, val6, val8, val9, val10, val11, val12, val13, val14, val15, val16, val17, val18, val19):
        self.w = val1 # Phonon frequencies
        self.v = val2 # Phonon eigenvectors
        self.atomsp = val4 # atom species (string)
        self.atomNs = val5 # atom species numbers (int)
        self.compound_name = val6 # compound_name (for plot saving)
        self.PBCTotal = val8
        self.PBCO = val9
        self.PBCA = val10
        self.PBCB = val11
        self.counter_TotalDOS = val12
        self.counter_partialDOS_O = val13
        self.counter_partialDOS_A = val14
        self.counter_partialDOS_B = val15
        self.w_TotalDOS = val16
        self.w_partialDOS_O = val17
        self.w_partialDOS_A = val18
        self.w_partialDOS_B = val19

    def plot_PBC(self, AB_first):
        freq = self.w.copy()
        #self.En /= np.sum(self.En) # Normalizing En
        #En = self.En.copy()
        TotalPBC = self.PBCTotal.copy()
        partialPBC = self.PBCO.copy()
        partialPBC_A = self.PBCA.copy()
        partialPBC_B = self.PBCB.copy()
        counter_TotalDOS = self.counter_TotalDOS.copy()
        counter_partialDOS_O = self.counter_partialDOS_O.copy()
        counter_partialDOS_A = self.counter_partialDOS_A.copy()
        counter_partialDOS_B = self.counter_partialDOS_B.copy()
        w_TotalDOS = self.w_TotalDOS.copy()
        w_partialDOS_O = self.w_partialDOS_O.copy()
        w_partialDOS_A = self.w_partialDOS_A.copy()
        w_partialDOS_B = self.w_partialDOS_B.copy()

        #modalEn_instance = En_hist_n_center(freq, En, 50)
        #[w_En, counter_En, w_En_center] = modalEn_instance.hist_for_En_n_center()
        #phononBC_instance = phononBC(50, freq, self.v, self.atomNs, self.atomsp)
        #[w_TotalDOS, counter_TotalDOS, TotalPBC] = phononBC_instance.TotalDOS() # Calculate Total DOS histogram
        #[w_partialDOS, counter_partialDOS_O, partialPBC] = phononBC_instance.partialDOS('O', AB_first) # Calculate partial DOS histogram
        #[w_partialDOS_A, counter_partialDOS_A, partialPBC_A] = phononBC_instance.partialDOS('A', AB_first)
        #[w_partialDOS_B, counter_partialDOS_B, partialPBC_B] = phononBC_instance.partialDOS('B', AB_first)
        counter_partialDOS_O /= np.sum(counter_TotalDOS) # Normalizing partial DOS first
        counter_TotalDOS /= np.sum(counter_TotalDOS) # Then, normalizing total DOS

        # plot series 1
        plt.figure(self.compound_name)
        plt.rcParams.update({'font.size': 13})
        #plt.fill_between(w_En, 0, self.smooth(counter_En,3), color="black", alpha=0.5, label='Phonon contributions')
        plt.fill_between(w_TotalDOS, 0, self.smooth(counter_TotalDOS,3), color="blue", alpha=0.5, label='Total phonon DOS')
        plt.fill_between(w_partialDOS_O, 0, self.smooth(counter_partialDOS_O,3), color="red", alpha=0.5, label='Partial pDOS')
        #plt.plot(freq,En, linestyle='none', marker='o', markerfacecolor='none', markersize=7, markeredgecolor='k', label='Phonon contributions')

        filename = "phonon_total_DOS_" + self.compound_name + ".txt"
        f_total_DOS = open(filename, "w")
        for i in np.arange(0,w_TotalDOS.shape[0],1):
            f_total_DOS.write("%f %f\n" %(w_TotalDOS[i], self.smooth(counter_TotalDOS,3)[i]))
        f_total_DOS.close()

        filename = "phonon_partial_DOS_" + self.compound_name + ".txt"
        f_partial_DOS = open(filename, "w")
        for i in np.arange(0,w_partialDOS_O.shape[0],1):
            f_partial_DOS.write("%f %f\n" %(w_partialDOS_O[i], self.smooth(counter_partialDOS_O,3)[i]))
        f_partial_DOS.close()

        minx = np.floor(np.min(freq))
        maxx = np.ceil(np.max(freq))
        miny = 0
        maxy = np.ceil(np.max(counter_TotalDOS)*20.)/20.

        #plt.plot([w_En_center,w_En_center],[miny,0.6*maxy],'k--')
        plt.plot([TotalPBC,TotalPBC],[miny,0.6*maxy],'b--')
        plt.plot([partialPBC,partialPBC],[miny,0.6*maxy],'r--')

        #plt.text(w_En_center, 0.63*maxy, 'En_center', fontsize=10, rotation = 'vertical', color = 'k')
        plt.text(TotalPBC, 0.63*maxy, 'tPBC', fontsize=10, rotation = 'vertical', color = 'b')
        plt.text(partialPBC, 0.63*maxy, 'pPBC', fontsize=10, rotation = 'vertical', color = 'r')
        
        plt.xlabel('Frequency (cm$^{-1}$)')
        plt.ylabel('Normlized excitation energy')
        plt.xlim(minx, maxx)
        plt.ylim(miny, maxy)

        plt.legend(fontsize=11)
        filename = 'phonon1_' + self.compound_name + '.png'
        plt.savefig(filename, bbox_inches = 'tight')
        #plt.show()
        plt.close()
        
    def smooth(self, y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

