void functest()
{
   auto f_combination = new TF1("f_combination",
    "[0] * (\
    (1-[4])*(3.*( [5]*2.*(1.-x**2) + [6]*(1.-2.*[3]*x+x**2) + (1.-[5]-[6])*(1.+2.*[3]*x+x**2) )/8.)\
    + [4]*(3.*( [1]*2.*(1.-x**2) + [2]*(1.-2.*[3]*x+x**2) + (1.-[1]-[2])*(1.+2.*[3]*x+x**2) )/8.))",
    -1,1);
   f_combination->SetParameters(10,0,0,0.16,1,0,0,0);
   f_combination->Draw();
   printf("%f\n",f_combination->Integral(-1,1));
}
