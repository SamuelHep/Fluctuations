

void TestCorCode()
{

    gSystem->Load("femtoLibrary/momentCode.so");
    FlucContainer * fluc = new FlucContainer();
    fluc->SetCumulants(1,1);
    fluc->SetCumulants(2,1);
    fluc->SetCumulants(3,1);
    fluc->SetCumulants(4,1);
    
    double m1 = fluc->GetMomentsFromCumulants_v2(1);
    double m2 = fluc->GetMomentsFromCumulants_v2(2);
    double m3 = fluc->GetMomentsFromCumulants_v2(3);
    double m4 = fluc->GetMomentsFromCumulants_v2(4);

    double c1 = fluc->GetCumulant(1);
    double c2 = fluc->GetCumulant(2);
    double c3 = fluc->GetCumulant(3);
    double c4 = fluc->GetCumulant(4);

    cout << "m1 = " << m1 << endl;
    cout << "m2 = " << m2 << endl;
    cout << "m3 = " << m3 << endl;
    cout << "m4 = " << m4 << endl;

    cout << "c1 = " << c1 << endl;
    cout << "c2 = " << c2 << endl;
    cout << "c3 = " << c3 << endl;
    cout << "c4 = " << c4 << endl;


}
