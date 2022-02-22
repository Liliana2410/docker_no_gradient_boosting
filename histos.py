from ROOT import TH1F


def histos():
	"""
	TH1F objetcs creation.
	- Convention of naming is: plot_kinematicOrTopologicalVariable_particle.
	- b1 an b2 are the leading b-jets and so on.
	Return: List of TH1F objetcs. 
	The last histo is plot_MET and the previous one is plot_PT_jets.
	"""
	#b-quark plots
	plot_PT_b1 = TH1F("Pt_b1", "Pt_b1", 100, 0.0, 1000.0)
	plot_ETA_b1 = TH1F("Eta_b1", "Eta_b1", 100, -5, 5) 
	plot_PHI_b1 = TH1F("Phi_b1", "Phi_b1", 100, -4, 4)  
	plot_PT_b2 = TH1F("Pt_b2", "Pt_b2", 100, 0.0, 1000.0) 
	plot_ETA_b2 = TH1F("Eta_b2", "Eta_b2", 100, -5, 5)
	plot_PHI_b2 = TH1F("Phi_b2", "Phi_b2", 100, -4, 4)

	#tau plots
	plot_PT_tau1 = TH1F("Pt_tau1", "Pt_tau1", 100, 0.0, 1000.0) 
	plot_ETA_tau1 = TH1F("Eta_tau1", "Eta_tau1", 100, -5, 5)
	plot_PHI_tau1 = TH1F("Phi_tau1", "Phi_tau1", 100, -4, 4)
	plot_PT_tau2 = TH1F("Pt_tau2", "Pt_tau2", 100, 0.0, 1000.0) 
	plot_ETA_tau2 = TH1F("Eta_tau2", "Eta_tau2", 100, -5, 5) 
	plot_PHI_tau2 = TH1F("Phi_tau2", "Phi_tau2", 100, -4, 4)

	#jet plots
	plot_PT_j1 = TH1F("Pt_j1", "Pt_j1", 100, 0.0, 1000.0) 
	plot_ETA_j1 = TH1F("Eta_j1", "Eta_j1", 100, -5, 5)
	plot_PHI_j1 = TH1F("Phi_j1", "Phi_j1", 100,-4, 4)
	plot_PT_j2 = TH1F("Pt_j2", "Pt_j2", 100, 0.0, 1000.0) 
	plot_ETA_j2 = TH1F("Eta_j2", "Eta_j2", 100, -5, 5)
	plot_PHI_j2 = TH1F("Phi_j2", "Phi_j2", 100,-4, 4)
	plot_PT_j3 = TH1F("Pt_j3", "Pt_j3", 100, 0.0, 1000.0) 
	plot_ETA_j3 = TH1F("Eta_j3", "Eta_j3", 100, -5, 5)
	plot_PHI_j3 = TH1F("Phi_j3", "Phi_j3", 100,-4, 4)
	plot_PT_j4 = TH1F("Pt_j4", "Pt_j4", 100, 0.0, 1000.0) 
	plot_ETA_j4 = TH1F("Eta_j4", "Eta_j4", 100, -5, 5)
	plot_PHI_j4 = TH1F("Phi_j4", "Phi_j4", 100,-4, 4)


	#Delta R
	plot_DeltaR_tau1tau2 = TH1F("DeltaR_tau1tau2", "DeltaR_tau1tau2", 100, 0, 10)
	plot_DeltaR_b1b2 = TH1F("DeltaR_b1b2", "DeltaR_b1b2", 100, 0, 10)
	plot_DeltaR_j1j2 = TH1F("DeltaR_j1j2", "DeltaR_j1j2", 100, 0, 10)
	plot_DeltaR_j1j3 = TH1F("DeltaR_j1j3", "DeltaR_j1j3", 100, 0, 10)
	plot_DeltaR_j1j4 = TH1F("DeltaR_j1j4", "DeltaR_j1j4", 100, 0, 10)
	plot_DeltaR_j2j3 = TH1F("DeltaR_j2j3", "DeltaR_j2j3", 100, 0, 10)
	plot_DeltaR_j2j4 = TH1F("DeltaR_j2j4", "DeltaR_j2j4", 100, 0, 10)
	plot_DeltaR_j3j4 = TH1F("DeltaR_j3j4", "DeltaR_j3j4", 100, 0, 10)

	#Delta Phi
	plot_DeltaPhi_tau1tau2 = TH1F("DeltaPhi_tau1tau2", "DeltaPhi_tau1tau2", 100, -4, 4)
	plot_DeltaPhi_b1b2 = TH1F("DeltaPhi_b1b2", "DeltaPhi_b1b2", 100, -4, 4)
	plot_DeltaPhi_j1j2 = TH1F("DeltaPhi_j1j2", "DeltaPhi_j1j2", 100, -4, 4)
	plot_DeltaPhi_j1j3 = TH1F("DeltaPhi_j1j3", "DeltaPhi_j1j3", 100, -4, 4)
	plot_DeltaPhi_j1j4 = TH1F("DeltaPhi_j1j4", "DeltaPhi_j1j4", 100, -4, 4)
	plot_DeltaPhi_j2j3 = TH1F("DeltaPhi_j2j3", "DeltaPhi_j2j3", 100, -4, 4)
	plot_DeltaPhi_j2j4 = TH1F("DeltaPhi_j2j4", "DeltaPhi_j2j4", 100, -4, 4)
	plot_DeltaPhi_j3j4 = TH1F("DeltaPhi_j3j4", "DeltaPhi_j3j4", 100, -4, 4)
	
	#Delta Pt
	plot_DeltaPt_tau1tau2 = TH1F("DeltaPt_tau1tau2", "DeltaPt_tau1tau2", 100, 0, 10)
	plot_DeltaPt_b1b2 = TH1F("DeltaPt_b1b2", "DeltaPt_b1b2", 100, 0, 10)
	plot_DeltaPt_j1j2 = TH1F("DeltaPt_j1j2", "DeltaPt_j1j2", 100, 0, 10)
	plot_DeltaPt_j1j3 = TH1F("DeltaPt_j1j3", "DeltaPt_j1j3", 100, 0, 10)
	plot_DeltaPt_j1j4 = TH1F("DeltaPt_j1j4", "DeltaPt_j1j4", 100, 0, 10)
	plot_DeltaPt_j2j3 = TH1F("DeltaPt_j2j3", "DeltaPt_j2j3", 100, 0, 10)
	plot_DeltaPt_j2j4 = TH1F("DeltaPt_j2j4", "DeltaPt_j2j4", 100, 0, 10)
	plot_DeltaPt_j3j4 = TH1F("DeltaPt_j3j4", "DeltaPt_j3j4", 100, 0, 10)

	#Reconstructed masses (j,j)
	plot_M_Zp = TH1F("M_Zp", "M_Zp", 100, 0.0, 3500.0)
	plot_M_j1j2 = TH1F("M_j1j2", "M_j1j2", 100, 0.0, 1000.0)
	plot_M_j1j3 = TH1F("M_j1j3", "M_j1j3", 100, 0.0, 1000.0)
	plot_M_j1j4 = TH1F("M_j1j4", "M_j1j4", 100, 0.0, 1000.0)
	plot_M_j2j3 = TH1F("M_j2j3", "M_j2j3", 100, 0.0, 1000.0)
	plot_M_j2j4 = TH1F("M_j2j4", "M_j2j4", 100, 0.0, 1000.0)
	plot_M_j3j4 = TH1F("M_j3j4", "M_j3j4", 100, 0.0, 1000.0)
	
	#Reconstructed masses (j,j,b1)
	plot_M_j1j2b1 = TH1F("M_j1j2b1", "M_j1j2b1", 100, 0.0, 1000.0)
	plot_M_j1j3b1 = TH1F("M_j1j3b1", "M_j1j3b1", 100, 0.0, 1000.0)
	plot_M_j1j4b1 = TH1F("M_j1j4b1", "M_j1j4b1", 100, 0.0, 1000.0)
	plot_M_j2j3b1 = TH1F("M_j2j3b1", "M_j2j3b1", 100, 0.0, 1000.0)
	plot_M_j2j4b1 = TH1F("M_j2j4b1", "M_j2j4b1", 100, 0.0, 1000.0)
	plot_M_j3j4b1 = TH1F("M_j3j4b1", "M_j3j4b1", 100, 0.0, 1000.0)
	
	#Reconstructed masses (j,j,b2)
	plot_M_j1j2b2 = TH1F("M_j1j2b2", "M_j1j2b2", 100, 0.0, 1000.0)
	plot_M_j1j3b2 = TH1F("M_j1j3b2", "M_j1j3b2", 100, 0.0, 1000.0)
	plot_M_j1j4b2 = TH1F("M_j1j4b2", "M_j1j4b2", 100, 0.0, 1000.0)
	plot_M_j2j3b2 = TH1F("M_j2j3b2", "M_j2j3b2", 100, 0.0, 1000.0)
	plot_M_j2j4b2 = TH1F("M_j2j4b2", "M_j2j4b2", 100, 0.0, 1000.0)
	plot_M_j3j4b2 = TH1F("M_j3j4b2", "M_j3j4b2", 100, 0.0, 1000.0)

	#StMet
	plot_StMet = TH1F("StMet", "StMet", 100, 0.0, 3000.0)

	#Pt of the most energetic taus
	plot_pt1 = TH1F("pt1", "pt1", 100, 0.0, 1000.0)
	plot_pt2 = TH1F("pt2", "pt2", 100, 0.0, 1000.0)

	#Fraction of pT/pL
	plot_PTvsPL = TH1F("PTvsPL", "PTvsPL", 100, 0.0, 1.0)

	#P of the two chosen taus
	plot_P_tau1 = TH1F("P_tau1", "P_tau1", 100, 0.0, 2000.0)
	plot_P_tau2 = TH1F("P_tau2", "P_tau2", 100, 0.0, 2000.0)

	#Aditional Elijah plots translated
	plot_M_tau1b1 = TH1F("M_tau1b1", "M_tau1b1", 100, 0.0, 1000.0)
	plot_M_tau1b2 = TH1F("M_tau1b2", "M_tau1b2", 100, 0.0, 1000.0)
	plot_M_b1b2 = TH1F("M_b1b2", "M_b1b2", 100, 0.0, 500.0)
	plot_M_tau2b1 = TH1F("M_tau2b1", "M_tau2b1", 100, 0.0, 1000.0)
	plot_M_tau2b2 = TH1F("M_tau2b2", "M_tau2b2", 100, 0.0 1000.0)
	plot_DeltaR_tau1b1 = TH1F("DeltaR_tau1b1", "DeltaR_tau1b1", 100, 0, 10)
	
	



	return (plot_PT_b1, plot_ETA_b1, plot_PHI_b1, plot_PT_b2, plot_ETA_b2, plot_PHI_b2, plot_PT_tau1, plot_ETA_tau1, plot_PHI_tau1, plot_PT_tau2, plot_ETA_tau2, plot_PHI_tau2, plot_PT_j1, plot_ETA_j1, plot_PHI_j1, plot_PT_j2, plot_ETA_j2, plot_PHI_j2, plot_PT_j3, plot_ETA_j3, plot_PHI_j3, plot_PT_j4, plot_ETA_j4, plot_PHI_j4, plot_DeltaR_tau1tau2, plot_DeltaR_b1b2, plot_DeltaR_j1j2, plot_DeltaR_j1j3, plot_DeltaR_j1j4, plot_DeltaR_j2j3, plot_DeltaR_j2j4, plot_DeltaR_j3j4, plot_DeltaPhi_tau1tau2, plot_DeltaPhi_b1b2, plot_DeltaPhi_j1j2, plot_DeltaPhi_j1j3, plot_DeltaPhi_j1j4, plot_DeltaPhi_j2j3, plot_DeltaPhi_j2j4, plot_DeltaPhi_j3j4, plot_DeltaPt_tau1tau2, plot_DeltaPt_b1b2, plot_DeltaPt_j1j2, plot_DeltaPt_j1j3, plot_DeltaPt_j1j4, plot_DeltaPt_j2j3, plot_DeltaPt_j2j4, plot_DeltaPt_j3j4, plot_M_Zp, plot_M_j1j2, plot_M_j1j3, plot_M_j1j4, plot_M_j2j3, plot_M_j2j4, plot_M_j3j4, plot_M_j1j2b1, plot_M_j1j3b1, plot_M_j1j4b1, plot_M_j2j3b1, plot_M_j2j4b1, plot_M_j3j4b1, plot_M_j1j2b2, plot_M_j1j3b2, plot_M_j1j4b2, plot_M_j2j3b2, plot_M_j2j4b2, plot_M_j3j4b2, plot_StMet, plot_pt1, plot_pt2, plot_PTvsPL, plot_P_tau1, plot_P_tau2, plot_M_tau1b1, plot_M_tau1b2, plot_M_b1b2, plot_M_tau2b1, plot_M_tau2b2, plot_DeltaR_tau1b1)

