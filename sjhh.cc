



int sjhh(){
	      //get the information form the first tree
              TFile *file = new TFile("events.root");
              TTree *tree = (TTree*)file->Get("tree");
             //get the information of no bec events                 
             // TFile *f3 = new TFile("nobec.root");
             // TTree *tree3=(TTree*)f3->Get("tree3");
             //creat a new tree to construct event blending
               TFile *new_file = new TFile("copy.root","recreate");
               TTree *tree2 = new TTree("tree2","copy of tree");
               Double_t  Q12,QR,eb,Qnobec;
               Double_t P1[4],P2[4],P3[4],PR[4],PW[4];
               Double_t mp  = 938.272013;
               Double_t mdelta=1.232;
               
               //TH1F *h1 =new TH1F("h1","Histogram 1",100,0,1);
               //TH1F *h2 =new TH1F("h2","Histogram 2",100,0,1);
                             

              Int_t bin_x=70;
               Int_t bin_y=70;
               Double_t min_x=0;
               Double_t max_x=1;
               Double_t min_y=1;
               Double_t max_y=3;
               TH2D *th_dalitz=new TH2D("th_dalitz","Dalitz",bin_x,min_x,max_x,bin_y,min_y,max_y);
               Double_t p12,p23,p13;










              
               tree->SetBranchAddress("P1",P1);
               tree->SetBranchAddress("P2",P2);
               tree->SetBranchAddress("eb",&eb);            
               tree->SetBranchAddress("P3",P3);
               tree->SetBranchAddress("Q12",&Q12);
              // tree3->SetBranchAddress("Qnobec",&Qnobec);
               
               tree2->Branch("PR",PR,"PR[4]/D");
               tree2->Branch("PW",PW,"PW[4]/D");
               tree2->Branch("P3",P3,"P3[4]/D");
               tree2->Branch("QR",&QR,"QR/D",100);
               tree2->Branch("Q12",&Q12);                              
               tree2->Branch("eb",&eb);        
             //  tree2->Branch("Qnobec",&Qnobec);      
               TLorentzVector pr,pw,MO;//pr ,pw is the P of two pio,and MO is p
               TLorentzVector target(0.0, 0.0, 0.0, mp*0.001);
               Double_t Er ;


               TRandom *randomGenerator = new TRandom();
               int a =randomGenerator->Integer(100)+1;


Int_t N=(Int_t)tree->GetEntries();


 for(Int_t i=0; i<N; i++){
     //tree3->GetEntry(i);
     tree->GetEntry(i);
    PR[0]=P1[0];  PR[1]=P1[1];  PR[2]=P1[2];  PR[3]=P1[3];
        //h1->Fill(Qnobec);
    MO.SetPxPyPzE(P3[1],P3[2],P3[3],P3[0]);
    tree->GetEntry((i+a) % tree->GetEntries());
    PW[0]=P2[0];  PW[1]=P2[1];  PW[2]=P2[2];  PW[3]=P2[3];
    TLorentzVector pr, pw,pl;// pr is the first pio,pw is the second ,pl is the weaker one
    pr.SetPxPyPzE(PR[1],PR[2],PR[3],PR[0]);
    pw.SetPxPyPzE(PW[1],PW[2],PW[3],PW[0]);   
    Er=eb;
      TLorentzVector beam(0.0, 0.0, Er,  Er );                
        TLorentzVector MC =beam+target-pr-pw;
  if (    pr.Mag()>pw.Mag())
      {pl=pw;}else{pl=pr;}
    
  TLorentzVector MT=MO + pl;//MT is the mass of two particl
                      
     //the following is the first condition, you don't know which one is the weaker ,so you need to have a test   
   if (abs(MC.Mag()-MO.Mag())<0.02){
        //the following is the second condition,you should remember to add the {} in the last,,
        // you don't know which one is the weaker ,so you need to have a test between pr and pw. 
               //  if(  abs(sqrt(pow(Pl[0],2)-pow(Pl[1],2)-pow(Pl[2],2)-pow(Pl[3],2))-mdelta)<0.05){ 
        TLorentzVector Psum = pr - pw;
        if(Psum.Mag2() > 0) { 
        QR = 0.0;
    } else { 
        QR = sqrt(-Psum.Mag2());
    }
                                                                                                                                                        
        p12=pow(PR[0]+PW[0],2)-pow(PR[1]+PW[1],2)-pow(PR[2]+PW[2],2)-pow(PR[3]+PW[3],2);
        p23=pow(PW[0]+P3[0],2)-pow(PW[1]+P3[1],2)-pow(PW[2]+P3[2],2)-pow(PW[3]+P3[3],2);
                    
                    th_dalitz->Fill(p12,p23);
        tree2->Fill();
        //h2->Fill(QR);
                //                                              }
	}
}

//the nest thing is to draw a graph ,you need to make it better,so you can have a good thing.
//I think you need to have a thing, you can ask gpt for help.It is not difficult,you can make it out very easily.
//to buy a bread and go to your home.

//TCanvas *c1 = new TCanvas("c1","Histograms",800,600);
     //      h1->Draw();
    //       h2->Draw("same");
           // 计算两个直方图的比值
  //TH1F *hratio = (TH1F*)h1->Clone("hratio");
  //hratio->Divide(h2);

  // 创建一个新画布并绘制比值直方图
  //TCanvas *c2 = new TCanvas("c2", "Ratio", 800, 600);
  //hratio->Draw();
  //h1->SaveAs("Qnobec.C");
  //h2->SaveAs("Qbec.C");
  //hratio->SaveAs("QC.C");

TCanvas *tc=new TCanvas("title","name",1920,1080);

tc->Divide(1,1);
tc->cd(1);
th_dalitz->Draw("zcol");
th_dalitz->GetXaxis()->SetTitle("M^{2}(#pi^{0},#pi^{0})");
th_dalitz->GetYaxis()->SetTitle("M^{2}(p,#pi^{0})");
th_dalitz->SetTitle(" 60%delta of Dalitz");
th_dalitz->SaveAs("Dalitz.C");


new_file->Write();
new_file->Close(); 	
return 0;

}
