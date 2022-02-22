import ROOT
from ROOT import TLorentzVector, TH1F
import histos
from collections import Counter
from math import sin, cos, sinh


def jet_reconstruction(jets, bs):

	global dijet_1
	global b_dijet_1
	global dijet_2
	global b_dijet_2
	global reconstructed_W1
	global reconstructed_W2
	#merged category of the first pair of hadronic jets coming from the top.
	global notMerged_1
	global partiallyMerged_1
	global fullyMerged_1
	#merged category of the second pair of hadronic jets coming from the top.
	global notMerged_2
	global partiallyMerged_2
	global fullyMerged_2
	#index of the first b-jet from the hadronic decaying top.
	global index_b1
	#index of the second b-jet from the hadronic decaying top.
	global index_b2
	#index of the two pairs of jets coming from the hadronic tops (pairs are always j1-j2 and j3-j4).
	global index_j1
	global index_j2
	global index_j3
	global index_j4

	#global Delta_PT_b_alpha_b_beta

	best_Err = 999999.9
	MW = 80.379
	Mt = 172.76

	index_j_temp1 = 0
	index_j_temp2 = 0
	index_j_temp3 = 0
	index_j_temp4 = 0
	
	index_b_temp1 = 0
	index_b_temp2 = 0
	
	#loop over the first group of particles (j_1, j_2, b_1).
	for i in range(len(jets)):
		for j in range(i + 1, len(jets)):
				for k in range(len(bs)):
				
					#Loop over the second group of particles (j_3, j_4, b_2).
					for n in range(len(jets)):
						if ((n != i) and (n != j)):
							for m in range(n + 1, len(jets)):
								if ((m != i) and (m != j)):
									for l in range(len(bs)):
										if (l != k):
							
											jtemp1 = jets[i]
											jtemp2 = jets[j]
											jtemp3 = jets[n]
											jtemp4 = jets[m]
											btemp1 = bs[k]
											btemp2 = bs[l]
												
											Err_1 = (abs((jtemp1 + jtemp2 + btemp1).M() - Mt))*100/Mt
											Err_2 = (abs((jtemp3 + jtemp4 + btemp2).M() - Mt))*100/Mt
											Err = Err_1 + Err_2
					
											#Selection criteria.
											if (Err < best_Err):
												best_Err = Err
												index_j_temp1 = i
												index_j_temp2 = j
												index_j_temp3 = n
												index_j_temp4 = m
												index_b_temp1 = k
												index_b_temp2 = l
					
	
	#Initialize the first group of particles.				
	index_j1 = index_j_temp1
	index_j2 = index_j_temp2	
	index_b1 = index_b_temp1
	
	dijet_1 = jets[index_j1] + jets[index_j2]		
	dr_dijet_1 = jets[index_j1].DeltaR(jets[index_j2])
	b_dijet_1 = jets[index_j1] + jets[index_j2] + bs[index_b1]
	dr_b_dijet_1 = dijet_1.DeltaR(bs[index_b1])
	
	#Initialize the second group of particles.
	index_j3 = index_j_temp3
	index_j4 = index_j_temp4	
	index_b2 = index_b_temp2
	
	dijet_2 = jets[index_j3] + jets[index_j4]		
	dr_dijet_2 = jets[index_j3].DeltaR(jets[index_j4])
	b_dijet_2 = jets[index_j3] + jets[index_j4] + bs[index_b2]
	dr_b_dijet_2 = dijet_2.DeltaR(bs[index_b2])
    		
	#Merged category for the first group of particles.
	if (dr_dijet_1 > 0.8):

	 	notMerged_1 = True
	 	partiallyMerged_1 = False
	 	fullyMerged_1 = False


   	else:

		notMerged_1 = False
		
        	if (dr_b_dijet_1 > 1.0):
        
			partiallyMerged_1 = True
			fullyMerged_1 = False
			reconstructed_W1 = dijet_1

		else:

			partiallyMerged_1 = False
			fullyMerged_1 = True
		
	#Merged category for the second group of particles.	
	if (dr_dijet_2 > 0.8):

	 	notMerged_2 = True
	 	partiallyMerged_2 = False
	 	fullyMerged_2 = False


   	else:

		notMerged_2 = False
		
        	if (dr_b_dijet_2 > 1.0):
        
			partiallyMerged_2 = True
			fullyMerged_2 = False
			reconstructed_W2 = dijet_2

		else:

			partiallyMerged_2 = False
			fullyMerged_2 = True

	return jets[index_j1], jets[index_j2], jets[index_j3], jets[index_j4], bs[index_b1], bs[index_b2]


def tau_reconstruction(taus):

	#index of the two taus (pairs of jets tagged as coming from taus), does NOT include the energy from the neutrinos (MET).
	global index_tau1
	global index_tau2

	best_dPt = 99999999.9
	for i in range(len(taus)):
		for j in range(i + 1, len(taus)):
			
			dPt = abs(taus[i].Pt() - taus[j].Pt())
			if (dPt < best_dPt):
			
				best_dPt = dPt
				index_tau1 = i
				index_tau2 = j
					
	return taus[index_tau1], taus[index_tau2]


def cross_cleaning(jets, bs, taus):

	particles = jets + bs
	temp = []
	for i in range(len(taus)):
		
		crossed = False
		for j in range(len(particles)):

			if (taus[i].DeltaR(particles[j]) < 0.3):
				crossed = True

		if (crossed == False):
			
			temp.append(taus[i])

	return temp
##############################	ESTO NO ############################
	
def frac_twobody_pTvspL(particle1, particle2):

	px1 = (particle1.Pt())*(cos(particle1.Phi()))
	py1 = (particle1.Pt())*(sin(particle1.Phi()))
	pz1 = (particle1.Pt())*(sinh(particle1.Eta()))
	
	px2 = (particle2.Pt())*(cos(particle2.Phi()))
	py2 = (particle2.Pt())*(sin(particle2.Phi()))
	pz2 = (particle2.Pt())*(sinh(particle2.Eta()))

	pT = ((px1 + px2)**2 + (py1 + py2)**2)**0.5
	pL = abs(pz1 + pz2)
	
	p = ((px1 + px2)**2 + (py1 + py2)**2 + (pz1 + pz2)**2)**0.5

	frac = pT/p

	return frac
####################################################################

def PT(TLV):
	return TLV.Pt()

def P(TLV):
	return TLV.P()
    
def histos_fill(plot, variable):
	return plot.Fill(variable)


def histos_Draw(plot):
	return plot.Draw('HISTOS')  


def histos_Write(plot):
	return plot.Write()


def histos_Reset(plot):
	return plot.Reset('ICESM') 


entries = ["Zprime_tata_350", "Zprime_tata_1000", "Zprime_tata_1500", "Zprime_tata_3000", "ttbarh", "ttbarZ"]
jobs = [10,10,10,10,10,10]

#Boolean value that determines if plots will be generated or not.
makeplots = True



if (makeplots == True):

	c1 = ROOT.TCanvas("c1", "Titulo")    # ROOT canvas

	#Creation of empty TH1F objects (empty ROOT histograms)
	plots = histos.histos()   # plots is a list of TH1F objects



arrs = []
print("Begin of Data Reading")



for n_signal, signal in enumerate(entries):

	arr1 = []

	#Counter used for saving in txt the information of the Gen-particles of the first 10 events. Only used when gentxt == True
	txt_counter = 0

	#List used to calculate the average total reconstructed PT/PL of the Z' (calculated with tau1 and tau2). Only used when checktaus == True.
	pt_p = []


	if (makeplots == True):
		f = ROOT.TFile(signal + ".root", "recreate")

	for ind in range(1, jobs[n_signal] + 1):

		directory = str("/disco4/SIMULACIONES/Liliana/fully_hadronic/" + signal + "/" + signal + "_" + str(ind) + "/Events/run_01/tag_1_delphes_events.root")
		File = ROOT.TChain("Delphes;1")
		File.Add(directory)
		Number = File.GetEntries()

		print("Signal: " + signal + "_" + str(ind))

		for i in range(Number):
			Entry = File.GetEntry(i)

			#Initializes particles lists.
			jets = []
			bs = []
			METs = []
			taus = []

			if (checktaus == True):
				tausbad = []

			EntryFromBranch_j = File.Jet.GetEntries()
			for j in range(EntryFromBranch_j):

				BTag = File.GetLeaf("Jet.BTag").GetValue(j)
				TauTag = File.GetLeaf("Jet.TauTag").GetValue(j)

				#searches for jets.
				if (BTag != 1 and TauTag != 1):
					jet = TLorentzVector()
					jet_PT, jet_Eta, jet_Phi, jet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
					jet.SetPtEtaPhiM(jet_PT, jet_Eta, jet_Phi, jet_M)
					jets.append(jet)

				#searches for b_jets.
				elif (BTag == 1 and TauTag != 1):
					bjet = TLorentzVector()
					bjet_PT, bjet_Eta, bjet_Phi, bjet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
					bjet.SetPtEtaPhiM(bjet_PT, bjet_Eta, bjet_Phi, bjet_M)
					bs.append(bjet)

				#searches for taus.
				elif (TauTag == 1 and BTag != 1):
					tau = TLorentzVector()
					tau_PT, tau_Eta, tau_Phi, tau_M, tau_Charge = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j),  File.GetLeaf("Jet.Charge").GetValue(j)
					tau.SetPtEtaPhiM(tau_PT, tau_Eta, tau_Phi, tau_M)
					taus.append(tau)

				if (checktaus == True):
					#searches for jets with both tags, bad taus (tauTag and bTag).
					if (TauTag == 1 and BTag == 1):
						tau = TLorentzVector()
						tau_PT, tau_Eta, tau_Phi, tau_M, tau_Charge = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j),  File.GetLeaf("Jet.Charge").GetValue(j)
						tau.SetPtEtaPhiM(tau_PT, tau_Eta, tau_Phi, tau_M)
						tausbad.append(tau)

			# MET (neutrinos).
			Total_MET = 0
			EntryFromBranch_MET = File.MissingET.GetEntries()
			for j in range(EntryFromBranch_MET):
				MET = TLorentzVector()
				MET_PT, MET_Eta, MET_Phi, MET_M  = File.GetLeaf("MissingET.MET").GetValue(j), File.GetLeaf("MissingET.Eta").GetValue(j), File.GetLeaf("MissingET.Phi").GetValue(j), 0.0
				MET.SetPtEtaPhiM(MET_PT, MET_Eta, MET_Phi, MET_M)
				METs.append(MET)
				Total_MET += MET_PT
      

			if (makeplots == True):
				
				plot_Cuts.AddBinContent(1)

				if (checktaus == True):

					#print(len(taus), len(tausbad))	

					plot_tausmalos.AddBinContent(1, len(taus))
					
					plot_tausmalos.AddBinContent(2, len(tausbad))

					a = plot_tausmalos.GetXaxis()
					a.SetBinLabel(1,"Real taus")
					a.SetBinLabel(2,"Fake taus")


      			#Checks if the event has the minimum theoretical expected number of particles (the count might be higher due to particle-detector effects).
			if (len(jets) >= 4 and len(bs) >= 2 and len(taus) >= 2):

				if (makeplots == True):
					plot_Cuts.AddBinContent(2)

				jets.sort(reverse = True, key=PT)     
				bs.sort(reverse = True, key=PT)
				taus.sort(reverse = True, key=PT)

				#Defines the real jets and bs from the event.
				j1, j2, j3, j4, b1, b2 = jet_reconstruction(jets, bs)
				realjets = [j1, j2, j3, j4]
				realbs = [b1, b2]

				#Makes the cross cleaning, only saving taus with dR not less than 0.3 to any other jet.
				realtaus = cross_cleaning(realjets, realbs, taus)

				#Checks if there are at least two real taus in the event (taus that passed the cross cleaning).
				if (len(realtaus) >= 2):

					if (makeplots == True):
						plot_Cuts.AddBinContent(3)

					tau1, tau2 = tau_reconstruction(realtaus)

			
	
					#~CONTINUING WITH THE EVENT.				
					#Pt requirement of the different jets (20 for taus is LHC minimum).
					if ((j1.Pt() > 30) and (j2.Pt() > 30) and (j3.Pt() > 30) and (j4.Pt() > 30) and (b1.Pt() > 30) and (j2.Pt() > 30) and (tau1.Pt() > 20) and (tau2.Pt() > 20)):

						if (makeplots == True):
							plot_Cuts.AddBinContent(4)

						#LHC required eta for taus.				
						if ((abs(tau1.Eta()) < 2.3) and (abs(tau2.Eta()) < 2.3)):

							if (makeplots == True):
								plot_Cuts.AddBinContent(5)

							#Minimum possible separation between both taus.
							if (tau1.DeltaR(tau2) > 0.3):

								if (training == True):
									row = np.array([tau1.Pt(), tau2.Pt(), tau1.DelptaPhi(tau2), tau1.DeltaR(tau2), (tau1 + tau2).M() + Total_MET, (j3 + j4).M(), j1.Pt(), j1.DeltaPhi(j2), j1.DeltaPhi(j3), j1.DeltaPhi(j4)])
            								arr1.append(row)

								#Checks if the plots should be filled and fills them.
								if (makeplots == True):

									plot_Cuts.AddBinContent(6)

									#Variables to be plotted.
									variables = [b1.Pt(), b1.Eta(), b1.Phi(), b2.Pt(), b2.Eta(), b2.Phi(), tau1.Pt(), tau1.Eta(), tau1.Phi(), tau2.Pt(), tau2.Eta(), tau2.Phi(), j1.Pt(), j1.Eta(), j1.Phi(), j2.Pt(), j2.Eta(), j2.Phi(), j3.Pt(), j3.Eta(), j3.Phi(), j4.Pt(), j4.Eta(), j4.Phi(), tau1.DeltaR(tau2), b1.DeltaR(b2), j1.DeltaR(j2), j1.DeltaR(j3), j1.DeltaR(j4), j2.DeltaR(j3), j2.DeltaR(j4), j3.DeltaR(j4), cos(tau1.DeltaPhi(tau2)), b1.DeltaPhi(b2), j1.DeltaPhi(j2), j1.DeltaPhi(j3), j1.DeltaPhi(j4), j2.DeltaPhi(j3), j2.DeltaPhi(j4), j3.DeltaPhi(j4) , abs(tau1.Pt() - tau2.Pt()), abs(b1.Pt() - b2.Pt()), abs(j1.Pt() - j2.Pt()), abs(j1.Pt() - j3.Pt()), abs(j1.Pt() - j4.Pt()), abs(j2.Pt() - j3.Pt()), abs(j2.Pt() - j4.Pt()), abs(j3.Pt() - j4.Pt()) , (tau1 + tau2).M() + Total_MET, (j1+j2).M(), (j1+j3).M(), (j1+j4).M(), (j2+j3).M(), (j2+j4).M(), (j3+j4).M(), (j1+j2+b1).M(), (j1+j3+b1).M(), (j1+j4+b1).M(), (j2+j3+b1).M(), (j2+j4+b1).M(), (j3+j4+b1).M(), (j1+j2+b2).M(), (j1+j3+b2).M(), (j1+j4+b2).M(), (j2+j3+b2).M(), (j2+j4+b2).M(), (j3+j4+b2).M(), j1.Pt() + j2.Pt() + j3.Pt() + j4.Pt() + b1.Pt() + b2.Pt() + tau1.Pt() + tau2.Pt() + Total_MET, taus[0].Pt(), taus[1].Pt(), frac_twobody_pTvspL(tau1, tau2), tau1.P(), tau2.P(), M(tau1 + b1), M(tau1 + b2), M(b1 + b2), M(tau2 + b1), M(tau2 + b2), tau1.DeltaR(b1)]

									b = plot_Cuts.GetXaxis()
									b.SetBinLabel(1,"Cut 0")
									b.SetBinLabel(2,"Cut 1")
									b.SetBinLabel(3,"Cut 2")
									b.SetBinLabel(4,"Cut 3")
									b.SetBinLabel(5,"Cut 4")
									b.SetBinLabel(6,"Cut 5")
									
									for i in range(len(plots)):
								  		histos_fill(plots[i], variables[i])
									
									if (checktaus == True):
										if (tau2.Pt() < 80.0):
											histos_fill(plot_PT_tausmenores, tau2.Pt())
											histos_fill(plot_ETA_tausmenores, tau2.Eta())
											histos_fill(plot_PHI_tausmenores, tau2.Phi())

											histos_fill(plot_DeltaR_tausmenores, tau1.DeltaR(tau2))
											histos_fill(plot_DeltaPhi_tausmenores, tau1.DeltaPhi(tau2))
										
										if (tau2.Pt() > 80.0):
											histos_fill(plot_PT_tausmayores, tau2.Pt())
											histos_fill(plot_ETA_tausmayores, tau2.Eta())
											histos_fill(plot_PHI_tausmayores, tau2.Phi())

											histos_fill(plot_DeltaR_tausmayores, tau1.DeltaR(tau2))
											histos_fill(plot_DeltaPhi_tausmayores, tau1.DeltaPhi(tau2))
	



	if (makeplots == True):

		#Drawing the histograms.
		plot_Cuts.Draw('HISTOS')

		for plot in plots:
			histos_Draw(plot)

		if (checktaus == True):
			plot_tausmalos.Draw('HISTOS')

			plot_PT_tausmenores.Draw('HISTOS')
			plot_PT_tausmayores.Draw('HISTOS')
			plot_ETA_tausmenores.Draw('HISTOS')
			plot_ETA_tausmayores.Draw('HISTOS')
			plot_PHI_tausmenores.Draw('HISTOS')
			plot_PHI_tausmayores.Draw('HISTOS')

			plot_DeltaR_tausmenores.Draw('HISTOS')
			plot_DeltaR_tausmayores.Draw('HISTOS')
			plot_DeltaPhi_tausmenores.Draw('HISTOS')
			plot_DeltaPhi_tausmayores.Draw('HISTOS')


		# Updating the canvas.
		c1.Update()


		# Writing the histograms.
		plot_Cuts.Write()

		for plot in plots:
			histos_Write(plot)

		if (checktaus == True):
			plot_tausmalos.Write()

			plot_PT_tausmenores.Write()
			plot_PT_tausmayores.Write()
			plot_ETA_tausmenores.Write()
			plot_ETA_tausmayores.Write()
			plot_PHI_tausmenores.Write()
			plot_PHI_tausmayores.Write()

			plot_DeltaR_tausmenores.Write()
			plot_DeltaR_tausmayores.Write()
			plot_DeltaPhi_tausmenores.Write()
			plot_DeltaPhi_tausmayores.Write()

		# Closing the ROOT file where the histos were saved.
		f.Close()

		  
		# Reseting the TH1F objects for its use in the next signal or bkg file.
		plot_Cuts.Reset('ICESM')

		for plot in plots:
			histos_Reset(plot)

		if (checktaus == True):
			plot_tausmalos.Reset('ICESM')

			plot_PT_tausmenores.Reset('ICESM')
			plot_PT_tausmayores.Reset('ICESM')
			plot_ETA_tausmenores.Reset('ICESM')
			plot_ETA_tausmayores.Reset('ICESM')
			plot_PHI_tausmenores.Reset('ICESM')
			plot_PHI_tausmayores.Reset('ICESM')

			plot_DeltaR_tausmenores.Reset('ICESM')
			plot_DeltaR_tausmayores.Reset('ICESM')
			plot_DeltaPhi_tausmenores.Reset('ICESM')
			plot_DeltaPhi_tausmayores.Reset('ICESM')


print("End of data reading")

print("Begin of Data Preparation")

#if (training == True):
