#include <iostream>
#include <vector>
#include <string>
#include <TH1F.h>
#include <fstream>
using namespace std;

void merge_histos()
{ 
  vector<std::string> names_files { "ttbarh.root","ttbarZ.root","Zprime_tata_350.root", "Zprime_tata_1000.root", "Zprime_tata_1500.root", "Zprime_tata_3000.root"};

  //vector<std::string> names_files_wo_ext {"$t\\bar{t}h$", "$t\\bar{t}t\\bar{t}$", "$t\\bar{t}b\\bar{b}/h$", "$WWZ$", "$m(Z')=350\\,GeV$", "$m(Z')=500\\,GeV$", "$m(Z')=750\\,GeV$", "$m(Z')=1000\\,GeV$", "$m(Z')=1500\\,GeV$", "$m(Z')=2000\\,GeV$", "$m(Z')=3000\\,GeV$", "$m(Z')=4000\\,GeV$"};

  vector<std::string> names { "ttbarh", "ttbarZ", "m(Z')= 350 GeV", "m(Z')= 1000 GeV", "m(Z')= 1500 GeV", "m(Z')= 3000 GeV"};    //for the legend names

  vector<std::string> plots {"PT_b1", "PT_b2", "PT_tau1", "PT_tau2"};

  vector<std::string> plots_names {"p_{T}(b_{1})", "p_{T}(b_{2})", "p_{T}(#tau_{1})","p_{T}(#tau_{2})"};

  vector<std::string> x_labels {"p_{T}(b_{1}) [GeV]", "p_{T}(b_{2}) [GeV]", "p_{T}(#tau_{1}) [GeV]", "p_{T}(#tau_{2}) [GeV]"};


  vector<int> colors {3, 7, 6, 5, 2, 4, 8, 9, 1, 43, 97, 38};
  //vector<int> colors {82, 5, 3, 7, 6, 5, 2, 4, 8, 9, 1, 43, 97, 38};   //linecolor

  vector<int> linestyles {1, 1, 10, 9, 8, 7, 6, 5, 4, 3};
  //vector<int> linestyles {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

  TList *l = new TList();

  // 
  for(int i=0; i<plots.size(); i++)
  {
    THStack *hs = new THStack("hs", plots_names[i].c_str());
    TCanvas *c2 = new TCanvas(plots[i].c_str(),"Histos",1280,1024);  
    Double_t x_1,x_2;

   
    x_1 = 0.65;
    x_2 = 0.85;
    

    // for legends: 
    auto legend = new TLegend(x_1,0.65,x_2,0.85);   // object of the class TLegend()

    for (int j=0; j<names.size(); j++)
    {    
      TFile f(names_files[j].c_str());
      TH1F *h = (TH1F*)f.Get(plots[i].c_str());
      h-> SetDirectory(0);
      h->SetLineColor(colors[j]);
      h->SetLineStyle(linestyles[j]);
      h->SetLineWidth(2);
      //h->GetMean();
      //h->GetStdDev();
      

      // to normalize the distribution:
      h->Scale(1.0/h->Integral());

      if((names_files[j] == "Zprime_tata_350.root") || (names_files[j] == "Zprime_tata_1500.root") || (names_files[j] == "Zprime_tata_1000.root") || (names_files[j] == "Zprime_tata_3000.root") )
      {
        legend->AddEntry(h,names[j].c_str(),"l");
      }

      else
      {
        legend->AddEntry(h,names[j].c_str(),"l");
      }
      
     
      legend->SetBorderSize(0);
      hs->Add(h);
    }   

    hs->Draw("NOSTACK HIST");
    hs->GetXaxis()->SetTitle(x_labels[i].c_str());
    hs->GetYaxis()->SetTitle("a.u.");
    legend->Draw();
    l->Add(c2);
    std::string filename ="imgs/"+ plots[i] + ".png";
    c2->SaveAs(filename.c_str());
  }

  TFile* Output = new TFile("joined.root", "RECREATE"); 
  l->Write();
  Output->Close();

}
