// --- March 11th, 2014 ---
//  1) firstly calculate correlation function, secondly nomarlized two Q ditributions to the same area
//  root -l ~/Desktop/1_researches/bec/ana/embec225/macro/show_cf_gen.C"(0,\"dt-phaseSpace.root\",\"dt-sigma.root\",0.0,0.7,\"figSigma\")"
int show_cf_gen(int plmode=0,  //0:Q and CF in different canvas;  1: in same canvas.
                //char fnobe[100]="nobec.root",
                //char fbe[100]="events.root",
                double fit1=0.0,
                double fit2=0.7
                //char figName[100]="fig16091601"
)
{
    
    //costant
    double fmtoGeV=5.07;   // 1fm = 5.07 GeV^-1
    //mode control
    //fitting control
    //variables
    Int_t  nPoints,iPoint=0;
    Double_t y1,y2,x,y,error;
    Long64_t  nentries=0;
    Int_t  nc,ic;
    char fileName[100];
    char srcFile[100];
    
    //cout << "bec file  :  " << fbe << endl;
    //cout << "nobec file:  " << fnobe << endl;
    cout << "fitting region: [" << fit1 << "," << fit2 << "]" << "GeV" << endl;
    
    // Q distribution free of BE effects
    Double_t Q12;
    //TFile *f1=new TFile(fnobe);
      //   TFile *f1=new TFile("nobec.root");
               TFile *f1=new TFile("copy.root");
// here to correct      TFile *f1=new TFile("copy.root");

    TH1D *hQnb=new TH1D("hQnb","hQnb", 40, 0, 1.0);
    //hQnb=(TH1D*) f1->Get("hQ");
  //  TTree *tnb=(TTree*) f1->Get("tree3");
             TTree *tnb=(TTree*) f1->Get("tree2");

    nentries = tnb->GetEntries();
         //cout << "nentries   " <<  nentries   <<endl;
  //  tnb->SetBranchAddress("Qnobec",&Q12);
  tnb->SetBranchAddress("QR",&Q12);

    for(Long64_t jentry=0; jentry<nentries; jentry++) {
        tnb->GetEntry(jentry);
        hQnb->Fill(Q12);
    }
    
    // Q distribution with BE effects
   // TFile *f2=new TFile(fbe);
       TFile *f2=new TFile("events.root");

    TH1D *hQb=new TH1D("hQb","hQb", 40, 0, 1.0);
    //hQb=(TH1D*) f2->Get("hQ");
    TTree *tb=(TTree*) f2->Get("tree");
    nentries = tb->GetEntries();
           tb->SetBranchAddress("Q12",&Q12);
    for(Long64_t jentry=0; jentry<nentries; jentry++) {
        tb->GetEntry(jentry);
        hQb->Fill(Q12);
    }
    
    
    // normalizing
    Double_t  nom1=hQb->Integral();
    Double_t  nom2=hQnb->Integral();
    //hQnb->Scale(nom1/nom2);
    
    iPoint=0;
    nPoints = hQb->GetNbinsX();
    TGraphErrors  *ratio2 = new TGraphErrors(nPoints);
    ratio2->SetTitle("Correlation Function");
    for (Int_t i=0;i<nPoints;i++) {
        y1=hQb->GetBinContent(i+1);
        y2=hQnb->GetBinContent(i+1);
        //cout << y1 << " " << y2 << endl;
        if (y1 >0.0 && y2 >0.0 && nom1>0.0 && nom2>0.0) {
            x=hQb->GetBinCenter(i+1);
            y= y1/y2*(nom2/nom1);
            // cout << "x,y" << x << " " << y << endl;
            error= y*sqrt(1./y1+1./y2);
            ratio2->SetPoint(iPoint, x , y );
            ratio2->SetPointError(iPoint, 0.0, error);
            iPoint ++;
        }
    }
    ratio2->Set(iPoint);
    ratio2->SetMarkerStyle(20);
    
    // ----------------- plot -------------------
    TPad *pad[10];
    TCanvas *canv[10];
    cout << "test" << endl;
    if(plmode==0) {
        ic=0;
        canv[ic] =new TCanvas("canv0","canv0",0,0,400,300);
        canv[ic]->cd(1); pad[0]=(TPad *) gPad;
        Int_t ic=1;
        canv[ic] =new TCanvas("canv1","canv1",500,0,400,300);
        canv[ic]->cd(1); pad[1]=(TPad *) gPad;
    }else if(plmode==1){
        Int_t ic=0;
        canv[ic] =new TCanvas("canv0","canv0",500,0,400,600);
        canv[ic]->Divide(1,2);
        canv[ic]->cd(1); pad[0]=(TPad *) gPad;
        canv[ic]->cd(2); pad[1]=(TPad *) gPad;
    }else{
        exit;
    }
    cout << "test2" << endl;
    // Q distribution
   hQnb->Scale(nom1/nom2);
    pad[0]->cd();
    hQb->SetTitle("Q (GeV)");
    hQb->SetLineWidth(3.0);
    hQb->GetXaxis()->SetTitle("Q (GeV)"); hQb->GetXaxis()->CenterTitle();hQb->GetXaxis()->SetLabelSize(0.06); hQb->GetXaxis()->SetTitleSize(0.07); // INPUT <<<<<<<<<<<<<<
    hQb->GetYaxis()->SetTitle("Counts"); hQb->GetYaxis()->CenterTitle();hQb->GetYaxis()->SetLabelSize(0.06); hQb->GetYaxis()->SetTitleSize(0.07); // INPUT <<<<<<<<<<<<<
    hQb->GetXaxis()-> SetNdivisions(5,5 ,0 , kTRUE);
    hQb->GetYaxis()-> SetNdivisions(5,5 ,0 , kTRUE);
    hQb->SetTitle("");
    hQb->SetMaximum(1.5*hQb->GetMaximum());
    hQb->SetStats(0);
    hQb->Draw("");
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.05);
    hQnb->SetLineColor(8);
    hQnb->SetLineWidth(3.0);
    hQnb->Draw("same");
    TLegend *legend = new TLegend(0.2,0.7,0.8,0.9);
    legend->AddEntry(hQb,"BEC");
    //legend->AddEntry(hQb,"#gamma p#rightarrow #Delta #pi^{0}#rightarrow#pi^{0}#pi^{0}p");
    //legend->AddEntry(hQb,"#gamma p#rightarrow #sigma p#rightarrow#pi^{0}#pi^{0}p");
    legend->AddEntry(hQnb,"MIX");
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.06);
    legend->Draw();
    
    // correlation function
    double Xmin=-0.01,Xmax=1.0;   // INPUT <<<<<<<<<<<<<<<
    double Ymin=0.0,Ymax=3.1;   // INPUT <<<<<<<<<<<<<<<
    TH1F *h=new TH1F("frame2","", 1000,Xmin,Xmax);
    h->SetStats(0); h->SetMinimum(Ymin); h->SetMaximum(Ymax);
    h->GetXaxis()->SetTitle("Q (GeV)"); h->GetXaxis()->CenterTitle();h->GetXaxis()->SetLabelSize(0.06); h->GetXaxis()->SetTitleSize(0.07);
    h->GetYaxis()->SetTitle("ratio"); h->GetYaxis()->CenterTitle();h->GetYaxis()->SetLabelSize(0.06); h->GetYaxis()->SetTitleSize(0.07);
    h->GetXaxis()-> SetNdivisions(10,5 ,0 , kTRUE);
    h->GetYaxis()-> SetNdivisions(5,5 ,0 , kTRUE);
    pad[1]->cd();  h->Draw();
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    // canv->cd(2); ratio2->Draw("p");
    pad[1]->cd(); ratio2->Draw("p");
    //Fitting
    TF1 *func;
    sprintf(fileName,"fit%3d",101);
    func = new TF1(fileName,"[0]*(1.0+[1]*TMath::Exp(-1*[2]*[2]*x*x))*(1.0)",0,0.95);
    func->SetParameter(0,0.95);
    func->SetParameter(1,0.5);
    func->SetParameter(2,2.0);
    func->SetParameter(3,1.9);
    ratio2->Fit(func,"","",fit1,fit2);
    TLatex *tex =new TLatex();
    tex->SetNDC(); tex->SetTextFont(62); tex->SetTextAlign(12);
    sprintf(fileName,"N= %2.2f#pm%2.2f", func->GetParameter(0), func->GetParError(0));
    tex->SetTextColor(2); tex->SetTextSize(0.06);
    tex->DrawLatex(0.63, 0.9, fileName);
     tex = new TLatex(Xmin+0.6*(Xmax-Xmin),0.9*Ymax,fileName);
    
    sprintf(fileName,"#lambda_{2} = %2.2f#pm%2.2f", func->GetParameter(1), func->GetParError(1));
    tex->SetTextColor(2); tex->SetTextSize(0.06);
    tex->DrawLatex(0.6, 2.15, fileName);
      tex = new TLatex(Xmin+0.6*(Xmax-Xmin),0.8*Ymax,fileName);
    
    sprintf(fileName,"r_{0}= %2.2f#pm%2.2f fm",
            fabs(func->GetParameter(2))/fmtoGeV, func->GetParError(2)/fmtoGeV);
    tex->SetTextColor(2); tex->SetTextSize(0.06);
    tex->DrawLatex(0.6, 2.5, fileName);
    tex = new TLatex(Xmin+0.6*(Xmax-Xmin),0.7*Ymax,fileName);
    
    sprintf(fileName,"#alpha= %2.2f#pm%2.2f GeV^{-1}",
            func->GetParameter(3), func->GetParError(3));
    tex->SetTextColor(2); tex->SetTextSize(0.06);
   // tex->DrawLatex(0.2, 0.72, fileName);
    
    tex->SetTextColor(2); tex->SetTextSize(0.09);
    //  tex->DrawLatex(0.2, 0.2, "(d)");
    pad[1]->Update();
    
    /*
     canv->cd(3);
     hQsm->SetLineWidth(3.0);
     hQsm->Draw();
     gPad->SetBottomMargin(0.15);
     gPad->SetLeftMargin(0.15);
     gPad->SetTopMargin(0.1);
     gPad->SetRightMargin(0.05);
     canv->Update();
     */
    // save canvas
    //input a text in this fig
    /*
     cout << "Please enter some words to identify this fig: "; cin >> fileName;
     tex->SetTextColor(2); tex->SetTextSize(0.07);
     canv->cd(1);  tex->DrawLatex(0.4, 0.96, fileName);
     canv->Update();
     */
    //cout << "File : " << figName << " .pdf/.eps/.gif/.C " << "have been created  " << endl;
    tex->SetTextColor(14); tex->SetTextSize(0.05);
    //pad[0]->cd(); tex->DrawLatex(0.05, 0.96, figName);
    if(plmode==0) {
        nc=2;
    }else if(plmode==1){
        nc=1;
    }else{
        exit;
    }
    for(Int_t ic=0;ic<nc; ic++){
        //sprintf(fileName,"%s_%02d.pdf",figName,ic);    canv[ic]->SaveAs(fileName);
       //sprintf(fileName,"%s_%02d.eps",figName,ic);    canv[ic]->SaveAs(fileName);
        //sprintf(fileName,"%s_%02d.gif",figName,ic);    canv[ic]->SaveAs(fileName);
        //sprintf(fileName,"%s_%02d.C",  figName,ic);    canv[ic]->SaveAs(fileName);
    }


//canv0->SaveAs("Q0.C");
//canv1->SaveAs("Q1.C");


    return 0;
}
