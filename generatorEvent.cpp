
Double_t fmtoGeV=5.07;  
Double_t mdelta=1.232;
Double_t wdelta=0.117;
Double_t msigma=0.500;
Double_t wsigma=0.500;
Double_t mf980=0.980;
Double_t wf980=0.040;
Double_t mp  = 938.272013;
Double_t mpi0= 134.9766;
Double_t mpi = 139.57018;
Double_t mgm = 0.0;




Double_t BECF(Double_t r,Double_t lambda,Double_t Q){
    return (1 + lambda*exp(-r*r*Q*Q));
}


Double_t BECF_q(Double_t r,Double_t lambda,Double_t q){
    return (1 + lambda*exp(-r*r*q*q));
}


Double_t BECF_qz_qr_p2(Double_t r,
                       Double_t lm,
                       Double_t qz, 
                       Double_t qr, 
                       Double_t p2)  
{
    return  (1.0+lm*exp(-qz*qz*r*r)*exp(-qr*qr*r*r/2.0)*ROOT::Math::cyl_bessel_j(0,qr/(p2*2.0)));
}


int generator(Int_t nepipip=30000, 
              Int_t nedelpi=30000, 
              Double_t eb1=1.15, 
              Double_t eb2=1.15, 
              Int_t ismear=0, 
              Int_t ics=0,  
              Int_t ibec=1,  
              Double_t becr = 0.8,  
              Double_t becst  = 0.5,  
              char *trName="events.root",
              char *dtClass="dtClass")
{
    gSystem->Load("libMathMore");
    printf ("%d pi0pi0p events will be generated ...\n", nepipip);
    printf ("%d Deltapi0 events will be generated ...\n", nedelpi);
    cout << "beam energy region: " << eb1 << " - " << eb2 << "" << " GeV" << endl;
    cout << "ismear:     " << ismear << endl;
    cout << "ics:       " << ics << endl;
    cout << "ibec:       " << ibec << endl;
    cout << "bec radius  :     " << becr << endl;
    cout << "bec strength:     " << becst << endl;
    cout << "tree name:  " << trName << endl;
    
    Int_t i,j,k;
    Double_t Q, weight;
    Int_t    BEC;
    Double_t eb;
    Double_t ranseed;
    Double_t Pi1[4], Pi2[4], Pi3[4], P1[4], P2[4], P3[4];
    Double_t theta1, theta2, theta3;
    Double_t phi1, phi2, phi3;
    Double_t M12, M23, M13, Q12, MM;
    Double_t Pm1[4],Pm2[4],Pm3[4],Pm4[4],Pm5[4];
   
    // check  mass of delta
    TH1D *hmdelta_free=new TH1D("hmdelta_free","",100,0,2.5);
    // check invariant mass of delta
    TH1D *hmdelta_ecut=new TH1D("hmdelta_ecut","",100,0,2.5);
    // check invariant mass of pi_0
    TH1D *hmdelta_picut=new TH1D("hmdelta_picut","",100,0,2.5);
    // check Q  relation function.
    TH1D *hQ=new TH1D("hQ","",100,0.0,1.0);
    // dalitz plot for delta 
    TH2D *hdalitz1=new TH2D("hdalitz1","",50,0,0.92, 50,0.98,3.03);
    // dalitz plot of delta ,saved only half events
    TH2D *hdalitz2=new TH2D("hdalitz2","",50,0.98,3.03, 50,0.98,3.03);
    // check two pi_0 invariant mass conbind 
    TH1D *hmppi=new TH1D("hmppi","",100,0.98,1.63);
    // check two pi_0 invariant mass M12 
    TH1D *hmppi1=new TH1D("hmppi1","",100,0.98,1.63);
    // check two pi_0 invariant mass M
    TH1D *hmppi2=new TH1D("hmppi2","",100,0.98,1.63);
    
    
    TFile f1(trName,"RECREATE");
    TTree *tree = new TTree( "tree", "phase space and bec events" );
    tree->Branch( "BEC", &BEC, "BEC/I"); 
    tree->Branch( "eb",  &eb, "eb/D"  );
    tree->Branch( "ran", &ranseed, "ran/D"  );
    tree->Branch( "P1",  P1, "P1[4]/D"  ); 
    tree->Branch( "P2",  P2, "P2[4]/D"  ); 
    tree->Branch( "P3",  P3, "P3[4]/D"  ); 
    tree->Branch("theta1", &theta1, "theta1/D",100 ); 
    tree->Branch("theta2", &theta2, "theta2/D",100 ); 
    tree->Branch("theta3", &theta3, "theta3/D",100 ); 
    tree->Branch("phi1",   &phi1,   "phi1/D",100 );
    tree->Branch("phi2",   &phi2,   "phi2/D",100 );
    tree->Branch("phi3",   &phi3,   "phi3/D",100 );
    tree->Branch("M12",    &M12,    "M12/D",100 );
    tree->Branch("M23",    &M23,    "M23/D",100 );
    tree->Branch("M13",    &M13,    "M13/D",100 );
    tree->Branch("Q12",    &Q12,    "Q12/D",100 );
    tree->Branch("MM",    &MM,    "MM/D",100 );
    tree->Branch( "Pm1", Pm1, "Pm1[4]/D"  ); 
    tree->Branch( "Pm2", Pm2, "Pm2[4]/D"  ); 
    tree->Branch( "Pm3", Pm3, "Pm3[4]/D"  ); 
    tree->Branch( "Pm4", Pm4, "Pm4[4]/D"  ); 
    tree->Branch( "Pm5", Pm5, "Pm5[4]/D"  ); 
    
    TLorentzVector target(0.0, 0.0, 0.0, mp*0.001);
    
    Double_t m[3];
    Double_t mf[2]={0.0, 0.0};
    TGenPhaseSpace event;
    TLorentzVector p1, p2, p3, q12, W, pdelta;
    TLorentzVector pm1, pm2,   pm3, pm4, pm5;
    TLorentzVector Psum;
    TRandom3 ran;
    Double_t r;
    Double_t BECweight;
    
    
    Int_t ievent=0;
    Double_t m3fs[3]={mpi0*0.001, mpi0*0.001, mp*0.001};
    while(ievent <nepipip) {
       
        Double_t Er = eb1+ran.Rndm()*(eb2-eb1);
        TLorentzVector beam(0.0, 0.0, Er,  Er );
        W = beam + target;
        Int_t ipermit=event.SetDecay(W,3,m3fs,"");
        if (ipermit !=1) continue;
        weight = event.Generate();
        p1=*(event.GetDecay(0));
        p2=*(event.GetDecay(1));
        p3=*(event.GetDecay(2));
        
        
        
        
        
        if(ibec==1) {
            q12 = p1-p2;
            if (q12.M2() <=0.0) {
                Q = - q12.M();
            }else{
                Q = 0.0;
            }
             BECweight = BECF(becr*fmtoGeV, becst, Q);
        }
        
        if(ibec==2) {
            TLorentzVector Ppi1, Ppi2;
            TLorentzVector Ppi1__, Ppi2__; 
            TLorentzVector Psum,Pppi1,Pppi2;
            TVector3 boost;
            Ppi1=p1; Ppi2=p2;
            Pppi1=p1+p3; Pppi2=p2+p3;
            if( fabs(Pppi1.M()-mdelta) < fabs(Pppi2.M()-mdelta) ) {
                boost.SetXYZ(-Pppi1.Px()/Pppi1.E(), -Pppi1.Py()/Pppi1.E(), -Pppi1.Pz()/Pppi1.E());
            } else {
                boost.SetXYZ(-Pppi2.Px()/Pppi2.E(), -Pppi2.Py()/Pppi2.E(), -Pppi2.Pz()/Pppi2.E());
            }
            Ppi1__=Ppi1; Ppi2__=Ppi2;
            Ppi1__.Boost(boost);
            Ppi2__.Boost(boost);
            Psum=Ppi1__-Ppi2__;
            BECweight = BECF_q(becr*fmtoGeV, becst, Psum.P());
        }
        
        
        if (ismear==1) {
            p1.SetRho(ran.Gaus(p1.Rho(),p1.Rho()*0.07));
            p1.SetTheta(ran.Gaus(p1.Theta(),0.09));
            p1.SetPhi(ran.Gaus(p1.Phi(),0.09));
            p1.SetE(sqrt(p1.Rho()*p1.Rho()+mpi0*mpi0*1.0e-6));
            p2.SetRho(ran.Gaus(p2.Rho(),p2.Rho()*0.07));
            p2.SetTheta(ran.Gaus(p2.Theta(),0.09));
            p2.SetPhi(ran.Gaus(p2.Phi(),0.09));
            p2.SetE(sqrt(p2.Rho()*p2.Rho()+mpi0*mpi0*1.0e-6));
            p3.SetRho(ran.Gaus(p3.Rho(),p3.Rho()*0.07));
            p3.SetTheta(ran.Gaus(p3.Theta(),0.09));
            p3.SetPhi(ran.Gaus(p3.Phi(),0.09));
            p3.SetE(sqrt(p3.Rho()*p3.Rho()+mp*mp*1.0e-6));
        }
        eb    = Er;
        ranseed   = ran.Rndm();
        P1[0] = p1.E();    P1[1] = p1.Px();    P1[2] = p1.Py();    P1[3] = p1.Pz();
        P2[0] = p2.E();    P2[1] = p2.Px();    P2[2] = p2.Py();    P2[3] = p2.Pz();
        P3[0] = p3.E();    P3[1] = p3.Px();    P3[2] = p3.Py();    P3[3] = p3.Pz();
        
        TLorentzVector Psum;
        theta1=p1.Theta(); phi1=p1.Phi();
        theta2=p2.Theta(); phi2=p2.Phi();
        theta3=p3.Theta(); phi3=p3.Phi();
        Psum = p1+p2; M12 = Psum.Mag();
        Psum = p2+p3; M23 = Psum.Mag();
        Psum = p1+p3; M13 = Psum.Mag();
        Psum = p1-p2; if(Psum.Mag2()> 0) { Q12=0.0;}else {Q12 = sqrt(-Psum.Mag2());}
        MM=p3.Mag();
        r=ran.Rndm();
       
        if (r <= weight){
            
            
            W = p1;
            event.SetDecay(W,2,mf,"");
            weight = event.Generate();
            pm1=*(event.GetDecay(0));
            pm2=*(event.GetDecay(1));
            W = p2;
            event.SetDecay(W,2,mf,"");
            weight = event.Generate();
            pm3=*(event.GetDecay(0));
            pm4=*(event.GetDecay(1));
            pm5=p3;
            Pm1[0]=pm1.E(); Pm1[1]=pm1.Px(); Pm1[2]=pm1.Py(); Pm1[3]=pm1.Pz();
            Pm2[0]=pm2.E(); Pm2[1]=pm2.Px(); Pm2[2]=pm2.Py(); Pm2[3]=pm2.Pz();
            Pm3[0]=pm3.E(); Pm3[1]=pm3.Px(); Pm3[2]=pm3.Py(); Pm3[3]=pm3.Pz();
            Pm4[0]=pm4.E(); Pm4[1]=pm4.Px(); Pm4[2]=pm4.Py(); Pm4[3]=pm4.Pz();
            Pm5[0]=pm5.E(); Pm5[1]=pm5.Px(); Pm5[2]=pm5.Py(); Pm5[3]=pm5.Pz();
            
            BEC = 0;
            if(ibec==0) { tree->Fill(); ievent++; hQ->Fill(Q12); continue;}
            
            r=ran.Rndm();
            if ( (BECweight/2.0) >r ){
                BEC =1;
                if(ibec!=0) { tree->Fill();  hQ->Fill(Q12); ievent++;}
            }
        }
    }
    
    
    
    ievent=0;
    while(ievent <nedelpi) {
        Double_t Er = eb1+ran.Rndm()*(eb2-eb1);
        
        m[1]=mpi0*0.001;
        TLorentzVector beam(0.0, 0.0, Er,  Er );
        W = beam + target;
        m[0]=ran.BreitWigner(mdelta,wdelta);
        hmdelta_free->Fill(m[0]);
        if(m[0]<0.0 || m[0]> W.M() ) continue;
       
        Int_t ipermit=event.SetDecay(W,2,m,"");
        if (ipermit !=1) continue;
        weight = event.Generate();
        pdelta=*(event.GetDecay(0));
        hmdelta_ecut->Fill(pdelta.M());
        p1=*(event.GetDecay(1));
        m[0]=mp*0.001;
        m[1]=mpi0*0.001;
       
        ipermit=event.SetDecay(pdelta,2,m,"");
        if (ipermit !=1) continue;
        hmdelta_picut->Fill(pdelta.M());
        weight = event.Generate();
        p3=*(event.GetDecay(0));
        p2=*(event.GetDecay(1));
        
        
        
        
        
        if(ibec==1) {
            q12 = p1-p2;
            if(q12.M2()<=0.0) { Q = - q12.M(); }else{Q = 0.0;}
            hQ->Fill(Q);
            BECweight = BECF(becr*fmtoGeV, becst, Q);
        }
        
        if(ibec==2 || ibec==3) {
            TLorentzVector Ppi1, Ppi2;
            TLorentzVector Ppi1__, Ppi2__, q_; 
            TLorentzVector Psum,Pppi1,Pppi2;
            TVector3 boost;
            TVector3 q_v3, Ppi1_v3, Ppi2_v3;
            Double_t q_z, q_r;
            Ppi1=p1; Ppi2=p2;
            Pppi1=p1+p3; Pppi2=p2+p3;
            boost.SetXYZ(-Pppi2.Px()/Pppi2.E(), -Pppi2.Py()/Pppi2.E(), -Pppi2.Pz()/Pppi2.E());
            Ppi1__=Ppi1; Ppi2__=Ppi2;
            Ppi1__.Boost(boost);
            Ppi2__.Boost(boost);
            q_=Ppi1__-Ppi2__;
            q_v3=q_.Vect();
            Ppi1_v3=Ppi1__.Vect();
            Ppi2_v3=Ppi2__.Vect();
            q_r=q_v3.Perp(Ppi2_v3);
            q_z=sqrt(q_.P()*q_.P()-q_r*q_r);
            if(ibec==2)  BECweight = BECF_q(becr*fmtoGeV, becst, q_.P());
            if(ibec==3)  BECweight = BECF_qz_qr_p2(becr*fmtoGeV, becst, q_z, q_r, Ppi2__.P());
            
        }
        
        if(ics==1){
            TLorentzVector Ppi1, Ppi2;
            TLorentzVector Ppi1__, Ppi2__; 
            TLorentzVector Psum,Pppi1,Pppi2;
            TVector3 boost;
            Ppi1=p1; Ppi2=p2;
            boost.SetXYZ(-pdelta.Px()/pdelta.E(), -pdelta.Py()/pdelta.E(), -pdelta.Pz()/pdelta.E());
            Ppi2__=Ppi2;
            Ppi2__.Boost(boost);
            Double_t CSweight = 2+3*sin(Ppi2__.Theta())*sin(Ppi2__.Theta());
            r=ran.Rndm();
            if ( (CSweight/5.0) < r ){
                continue;
            }
        }
        
        
        
        
        if (ismear==1) {
            p1.SetRho(ran.Gaus(p1.Rho(),p1.Rho()*0.07));
            p1.SetTheta(ran.Gaus(p1.Theta(),0.09));
            p1.SetPhi(ran.Gaus(p1.Phi(),0.09));
            p1.SetE(sqrt(p1.Rho()*p1.Rho()+mpi0*mpi0*1.0e-6));
            p2.SetRho(ran.Gaus(p2.Rho(),p2.Rho()*0.07));
            p2.SetTheta(ran.Gaus(p2.Theta(),0.09));
            p2.SetPhi(ran.Gaus(p2.Phi(),0.09));
            p2.SetE(sqrt(p2.Rho()*p2.Rho()+mpi0*mpi0*1.0e-6));
            p3.SetRho(ran.Gaus(p3.Rho(),p3.Rho()*0.07));
            p3.SetTheta(ran.Gaus(p3.Theta(),0.09));
            p3.SetPhi(ran.Gaus(p3.Phi(),0.09));
            p3.SetE(sqrt(p3.Rho()*p3.Rho()+mp*mp*1.0e-6));
        }
        TLorentzVector Psum;
        eb    = Er;
        ran   = ran.Rndm();
        P1[0] = p1.E();    P1[1] = p1.Px();    P1[2] = p1.Py();    P1[3] = p1.Pz();
        P2[0] = p2.E();    P2[1] = p2.Px();    P2[2] = p2.Py();    P2[3] = p2.Pz();
        P3[0] = p3.E();    P3[1] = p3.Px();    P3[2] = p3.Py();    P3[3] = p3.Pz();
        theta1=p1.Theta(); phi1=p1.Phi();
        theta2=p2.Theta(); phi2=p2.Phi();
        theta3=p3.Theta(); phi3=p3.Phi();
        Psum = p1+p2; M12 = Psum.Mag();
        Psum = p2+p3; M23 = Psum.Mag();
        Psum = p1+p3; M13 = Psum.Mag();
        Psum = p1-p2; if(Psum.Mag2()> 0) { Q12=0.0;}else {Q12 = sqrt(-Psum.Mag2());}
        r=ran.Rndm();
        MM=p3.M();
        if (r <= weight){
            
            W = p1;
            event.SetDecay(W,2,mf,"");
            weight = event.Generate();
            pm1=*(event.GetDecay(0));
            pm2=*(event.GetDecay(1));
            W = p2;
            event.SetDecay(W,2,mf,"");
            weight = event.Generate();
            pm3=*(event.GetDecay(0));
            pm4=*(event.GetDecay(1));
            pm5=p3;
            Pm1[0]=pm1.E(); Pm1[1]=pm1.Px(); Pm1[2]=pm1.Py(); Pm1[3]=pm1.Pz();
            Pm2[0]=pm2.E(); Pm2[1]=pm2.Px(); Pm2[2]=pm2.Py(); Pm2[3]=pm2.Pz();
            Pm3[0]=pm3.E(); Pm3[1]=pm3.Px(); Pm3[2]=pm3.Py(); Pm3[3]=pm3.Pz();
            Pm4[0]=pm4.E(); Pm4[1]=pm4.Px(); Pm4[2]=pm4.Py(); Pm4[3]=pm4.Pz();
            Pm5[0]=pm5.E(); Pm5[1]=pm5.Px(); Pm5[2]=pm5.Py(); Pm5[3]=pm5.Pz();
            BEC = 0;
            if(ibec==0) { tree->Fill(); ievent++; continue;}
            
            r=ran.Rndm();
            if ( (BECweight/2.0) >r ){
                BEC =1;
                if(ibec!=0) { tree->Fill(); ievent++;}
            }
            
            hdalitz1->Fill(M12*M12,M13*M13);
            hdalitz1->Fill(M12*M12,M23*M23);
            if(ievent%2==0) {
                hdalitz2->Fill(M13*M13,M23*M23);
            }else{
                hdalitz2->Fill(M23*M23,M13*M13);
            }
            hmppi->Fill(M13);
            hmppi->Fill(M23);
            hmppi1->Fill(M23);
            hmppi2->Fill(M13);
        }
    }
    f1.cd();
    hmdelta_free->GetXaxis()->SetTitle("Mass of #Delta (MeV)");   
    hmdelta_free->GetXaxis()->CenterTitle(true);
    hmdelta_free->GetYaxis()->SetTitle("Counts");
    hmdelta_free->GetYaxis()->CenterTitle();
    hmdelta_free->SetLineColor(kBlue);
    hmdelta_free->SetTitle("Distribution of #Delta Mass");
    hmdelta_free->SetLineWidth(2);
    hmdelta_free->Write();
    hmdelta_free->Draw();



    hmdelta_ecut->GetXaxis()->SetTitle("Invariant Mass of #Delta (MeV)");
    hmdelta_ecut->GetXaxis()->CenterTitle();
    hmdelta_ecut->GetYaxis()->SetTitle("Counts");
    hmdelta_ecut->GetYaxis()->CenterTitle();
    hmdelta_ecut->SetLineColor(kBlue);
    hmdelta_ecut->SetTitle("Distribution of #Delta Invariant Mass");
    hmdelta_ecut->SetLineWidth(2);
    hmdelta_ecut->SetFillColor(kYellow);
    hmdelta_ecut->Write();
    hmdelta_ecut->Draw();


    hmdelta_picut->GetXaxis()->SetTitle("Invariant Mass of #Delta (MeV)");
    hmdelta_picut->GetXaxis()->CenterTitle();
    hmdelta_picut->GetYaxis()->SetTitle("Counts");
    hmdelta_picut->GetYaxis()->CenterTitle();
    hmdelta_picut->SetLineColor(kBlue);
    hmdelta_picut->SetTitle("Distribution of #Delta Invariant Mass");
    hmdelta_picut->SetLineWidth(2);
    hmdelta_picut->SetFillColor(kYellow);
    hmdelta_picut->Write();
    hmdelta_picut->Draw();
    


    hQ->Write();
    hmppi->Write();
    hmppi1->Write();
    hmppi2->Write();
    hdalitz1->Write();
    hdalitz2->Write();
    tree->Write();
    tree->MakeClass(dtClass);
    f1.Close();  
     
     
     
     
     
    
       
       
       
    
       
       
        
    
     
     
     
    
  
    
      
    
    return 0;
}
