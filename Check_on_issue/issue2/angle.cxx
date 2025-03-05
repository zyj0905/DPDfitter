void print_angles1(TLorentzVector T4_1,TLorentzVector T4_2,TLorentzVector T4_3){
    double alpha = (-T4_1).Phi();
    double beta = (-T4_1).Theta();
    cout<<"The correct code, "<<"beta: "<<beta<<", alpha:"<<alpha<<endl;
    T4_1.RotateZ(-alpha); T4_2.RotateZ(-alpha); T4_3.RotateZ(-alpha);
    T4_1.RotateY(-beta); T4_2.RotateY(-beta); T4_3.RotateY(-beta);

    TVector3 boostVector = -(T4_2+T4_3).BoostVector();
    T4_1.Boost(boostVector); T4_2.Boost(boostVector); T4_3.Boost(boostVector);
    double theta23 = T4_2.Theta();
    double phi23 = T4_2.Phi();
    cout<<"The correct code, "<<"theta23:  "<<theta23<<", phi23:"<<phi23<<endl;
}

void print_angles2(TLorentzVector T4_1,TLorentzVector T4_2,TLorentzVector T4_3){
    double alpha = (-T4_1).Phi();
    double beta = (-T4_1).Theta();
    cout<<"My previous code, "<<"beta: "<<beta<<", alpha:"<<alpha<<endl;

    TVector3 plane23 = T4_3.Vect().Cross(T4_2.Vect());
    TVector3 Z_axis(0,0,1);
    TVector3 plane1Z = Z_axis.Cross(T4_1.Vect());
    double phi23 = ((plane1Z.Cross(plane23)).Dot(T4_1.Vect())>0? 1.0 : -1.0) * (plane1Z).Angle(plane23);

    TVector3 boostVector = -(T4_2+T4_3).BoostVector();
    T4_1.Boost(boostVector); T4_2.Boost(boostVector); T4_3.Boost(boostVector);
    double theta23 = (-T4_1.Vect()).Angle(T4_2.Vect());
    cout<<"My previous code, "<<"theta23:  "<<theta23<<", phi23:"<<phi23<<endl;
}

void angle(){
    TLorentzVector Tp4_Dst,Tp4_D,Tp4_Pi;
    TLorentzVector T4_1,T4_2,T4_3;
    Tp4_Dst.SetPxPyPzE(0.0570074,-0.026685,0.0479813,2.0085299);
    Tp4_D.SetPxPyPzE(-0.108349,-0.056907,-0.295222,1.8967276);
    Tp4_Pi.SetPxPyPzE(0.1034199, 0.0873037,0.2482869,0.3153471);

    TVector3 boostVector = -(Tp4_Dst+Tp4_D+Tp4_Pi).BoostVector();
    Tp4_Dst.Boost(boostVector);
    Tp4_D.Boost(boostVector);
    Tp4_Pi.Boost(boostVector);

    //Tp4_Dst.Print();
    //Tp4_D.Print();
    //Tp4_Pi.Print();

    //Here assign different reference chain

    //Dst pi D
    cout<<"Dst(Dpi) as reference chain:"<<endl;
    T4_1 = Tp4_Dst; T4_2 = Tp4_Pi; T4_3 = Tp4_D;
    print_angles1(T4_1,T4_2,T4_3);
    print_angles2(T4_1,T4_2,T4_3);

    //Pi D Dst
    cout<<"Pi(DDst) as reference chain:"<<endl;
    T4_1 = Tp4_Pi; T4_2 = Tp4_D; T4_3 = Tp4_Dst;
    print_angles1(T4_1,T4_2,T4_3);
    print_angles2(T4_1,T4_2,T4_3);

    //D Dst pi
    cout<<"D(Dstpi) as reference chain:"<<endl;
    T4_1 = Tp4_D; T4_2 = Tp4_Dst; T4_3 = Tp4_Pi;
    print_angles1(T4_1,T4_2,T4_3);
    print_angles2(T4_1,T4_2,T4_3);


}
