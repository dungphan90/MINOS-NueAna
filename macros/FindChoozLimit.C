float FindChoozLimit(float dm2){

 std::ifstream ins;
 ins.open("ChoozPoints.txt");
 float d1, d2, d3;
 
 float chooz_deltam[200];
 float chooz_sinsq[200];
 
 ins>>d1>>d2>>d3;
 
 int count = 0;
 while(!ins.eof()){
   chooz_deltam[count] = d2;
   chooz_sinsq[count] = d1;
   count++;
   ins>>d1>>d2>>d3;
 }
 
  for(int i=1;i<200;i++){
     float choozdm = chooz_deltam[i];
     float choozdmprev= chooz_deltam[i-1];
     if(dm2>choozdmprev&&dm2<choozdm){
	float m= (chooz_sinsq[i]-chooz_sinsq[i-1])/(choozdm-choozdmprev);
        cout<<choozdm<<"   "<<choozdmprev<<endl;
        cout<<m<<"  "<<chooz_sinsq[i-1]<<"  "<<(dm2-choozdmprev)<<"  "<<dm2<<endl;
        
        if(m == 0) return chooz_sinsq[i];
	float choozsue3 = chooz_sinsq[i-1]+(dm2-choozdmprev)*m;
	return choozsue3;
     }
  }
  return 0;
}
