clc
clear

sigma=10e-9; 
lambda=1; 
L=18*sigma;  


R= 4.0*sigma;  

nuL= (7.68*10^-3)/sigma^2; 
nu= nuL/L; 



%% N1 % Equation 9 in the paper
lambda_t= 2*pi*nuL;
f_lambda=  ((sqrt(lambda^2+1)-lambda)*lambda)+(log((sqrt(lambda^2+1))+lambda));
N1=lambda_t* f_lambda *  R^2;


%% N2  % Equation A2
lambda_e= (4/3)*pi*nu;
g_lambda= 1.5*lambda;
N2= lambda_e* g_lambda* R^3;
%% N4a % Equations are given in Appendix A
m= 9;
points= m+1;
f_1=[0]; 
f_2=[0];  
f_3=[0]; 
 

counter_r=1;
 
for r = linspace(R, R*sqrt(1+lambda^2),points) 

    alpha1=atan(1/lambda);    
    theta=asin (R/r); 
    alpha2_min=theta-alpha1;
    
    alpha2max_cond= sqrt (((1+lambda^2)*R^2)-(r^2)); 
    if ((2*L)>= alpha2max_cond )
        alpha2_max=acos(r/(sqrt(1+lambda^2)*R)); 
    else
        alpha2_max=asin((R*r+(2*L*sqrt((4*L^2)+(r^2-R^2))))/(r^2+(4*L^2)))-alpha1;
    end



        counter_alpha=1;
        for alpha2= linspace(alpha2_min, alpha2_max,points) 
            AC=sqrt(abs(((1+lambda^2)*R^2)+(r^2)-(2*R*r*cos(alpha2)*sqrt(1+lambda^2)))); 
            if ((2*L)>=AC )
               eta_max=atan(((sqrt(1+lambda^2)*R)-(r*cos(alpha2)))/(r*sin(alpha2)))- alpha2; 
            else
                
               eta_max= acos(((r*sin(alpha1+alpha2))-R)/(2*L))-(alpha1+alpha2) ;
            end

                if (isnan(eta_max)) 
                    eta_max=0;
                end
            
            counter_eta=1;
            for eta=linspace(0, eta_max,points)   
                
                AP=sqrt(([(r^2*sin(alpha1+alpha2)^2)-R^2]*[(r^2*sin(alpha1+alpha2)^2)-(R^2*sin(alpha1+alpha2+eta)^2)])/(r^2*[sin(alpha1+alpha2)^2]*[cos(alpha1+alpha2+eta)^2]));
                if ((2*L)>=AP)
                   beta_max= acot((sqrt((r^2*sin(alpha1+alpha2)^2)-R^2))/(R*cos(alpha1+alpha2+eta))); 
                else
                   beta_max= acos (((-r*sin(alpha1+alpha2)*cos(alpha1+alpha2+eta))+sqrt((r^2*sin(alpha1+alpha2)^2)+(sin(alpha1+alpha2+eta)^2*(4*L^2-R^2))))/(2*L*sin(alpha1+alpha2+eta)^2)); 
                end
                f_3(counter_eta)=sin(beta_max)/pi; 
                counter_eta=counter_eta+1;
            end
            int3= (eta_max/(2*m))*(f_3(1)+sum(2*f_3(2:end-1))+f_3(end)); 
            f_2(counter_alpha)=sin(alpha1+alpha2)*int3;
            counter_alpha=counter_alpha+1;
            
        end
        
        int2= ((alpha2_max-alpha2_min)/(2*m))* (f_2(1)+sum(2*f_2(2:end-1))+f_2(end)); 
        r;
        f_1(counter_r)=3*r^2*int2;
        counter_r=counter_r+1;
end
N4a= ((R*sqrt(1+lambda^2)-R)/(2*m))* (f_1(1)+sum(2*f_1(2:end-1))+f_1(end));
N4a= real (N4a*lambda_e);


%% N4b
m= 19;
points= m+1;
f_1=[0]; 
f_2=[0]; 
f_3=[0]; 
 

counter_r=1;

for r = linspace(lambda*R, R*sqrt(1+lambda^2),points)
    lambda_max= acos((lambda*R)/r);

    lambda_c= asin ([r-sqrt(r^2-[(1+lambda^2)*(r^2-[lambda*R]^2)])]/[R*(1+lambda^2)]);
    CE=[R-r*sin(lambda_c)]/cos(lambda_c);

    if ((2*L)>= CE)
      lambda_min= asin ([r-sqrt(r^2-[(1+lambda^2)*(r^2-lambda^2*R^2)])]/[R*(1+lambda^2)])  ;
    else
       lambda_min=asin ([(r*sqrt(4*L^2+r^2-[lambda*R]^2))-(2*L*lambda*R)]/[4*L^2+r^2]) ; 
    end
    counter_lada=1;  
    for lada= linspace(lambda_min, lambda_max,points) 
     CD= sqrt((r*cos(lada)-lambda*R)^2+(R-r*sin(lada))^2);
     
     if ((2*L)>=CD)
       eta_max=lada-atan (((r*cos(lada))-(lambda*R))/(R-r*sin(lada))) ; 
     else
       eta_max=lada-asin ((r*cos(lada)-lambda*R)/(2*L));
     end
     eta=linspace(0, eta_max,points);
     int3=trapz(eta,sin(eta)/pi);
     
      f_2(counter_lada)=sin(lada)*int3;
      counter_lada=counter_lada+1;
    end

    int2= ((lambda_max-lambda_min)/(2*m))* (f_2(1)+sum(2*f_2(2:end-1))+f_2(end));  
    f_1(counter_r)=3*r^2*int2;
    counter_r=counter_r+1;
end

N4b= ((R*(sqrt(1+lambda^2)-lambda))/(2*m))* (f_1(1)+sum(2*f_1(2:end-1))+f_1(end));
N4b= real (N4b*lambda_e);




%% N3a
M= 9;   
points= M+1;
f_1=[0]; 
f_2=[0]; 
f_3=[0]; 



counter_r=1;
for r = linspace(R, R*sqrt(1+lambda^2),points) 

    alpha1=atan(1/lambda);    
    theta=asin (R/r); 
    alpha2_min=theta-alpha1; 
    
    alpha2max_cond= sqrt (((1+lambda^2)*R^2)-(r^2));
    if ((2*L)>= alpha2max_cond )
        alpha2_max=acos(r/(sqrt(1+lambda^2)*R));
    else
        alpha2_max=asin((R*r+(2*L*sqrt((4*L^2)+(r^2-R^2))))/(r^2+(4*L^2)))-alpha1;
    end
    counter_alpha=1;
    
    
        for alpha2= linspace(alpha2_min, alpha2_max,points)  
            mmin=  (r*sin(alpha1+alpha2)-R)/(cos (alpha1+alpha2));  
  
            mmin= abs (mmin);  
            counter_m=1;
        
            
            for m=linspace(mmin, 2*L,points)  
                alpha_c=acos((r^2-R^2)/(lambda*R*r))-alpha1; 
                
                if alpha2 <= alpha_c
                
                    AP= sqrt (abs ([(r^2*sin(alpha1+alpha2)^2-R^2)*(r^2-R^2)]/[r*cos(alpha1+alpha2)]^2) ); 
             
                    if AP<=m
                     beta= acot ((sqrt(r^2*sin(alpha1+alpha2)^2-R^2))/[R*cos(alpha1+alpha2)])  ; 
                    else
                      
                     beta= acos (([-r*cos(alpha1+alpha2)]+sqrt(r^2+(m^2-R^2)))/(m*sin(alpha1+alpha2))); 
                     
                         if ([-r*cos(alpha1+alpha2)]+sqrt(r^2+(m^2-R^2)))==0 && (m*sin(alpha1+alpha2))==0 
                         beta= acos (0);    
                         end
                    end
                    
                else
                    AQ= sqrt((R^2*[1+lambda^2])-r^2);
                    if AQ>m
                     beta= acos (([-r*cos(alpha1+alpha2)]+sqrt(r^2+(m^2-R^2)))/(m*sin(alpha1+alpha2)));
                    
                    else
                     beta= atan ((sqrt(R^2-r^2-[(1+lambda^2)*(R*cos(alpha1+alpha2))^2]+[2*lambda*R*r*cos(alpha1+alpha2)]))/[lambda*R-r*cos(alpha1+alpha2)]);    
                    end    
                end
          
                f_3(counter_m)=beta / (pi*L); 
                counter_m=counter_m+1;

      
            end 

            int3= ((2*L-mmin)/(2*M))*(f_3(1)+sum(2*f_3(2:end-1))+f_3(end)); 

            f_2(counter_alpha)=sin(alpha1+alpha2)*int3;
            counter_alpha=counter_alpha+1;
        end

        int2= ((alpha2_max-alpha2_min)/(2*M))* (f_2(1)+sum(2*f_2(2:end-1))+f_2(end)); 
        f_1(counter_r)=2*r*int2;
        counter_r=counter_r+1;
end

N3a= ((R*sqrt(1+lambda^2)-R)/(2*M))* real(f_1(1)+sum(2*f_1(2:end-1))+f_1(end));
N3a= real (N3a*lambda_t);




%% N3b

m= 9;
points= m+1;
f_1=[0]; 
f_2=[0];  
f_3=[0]; 
counter_r=1;
epsilon= ((R*sqrt(1+lambda^2))-(lambda*R))/10;

for r = linspace(lambda*R+epsilon, epsilon+ R*sqrt(1+lambda^2),points)  

   lada_max= acos ((lambda*R)/r);  
   lada_c= asin ([r-sqrt(r^2-[(1+lambda^2)*(r^2-lambda^2*R^2)])]/[R*(1+lambda^2)]); 
   CE=(R-r*sin(lada_c))/cos(lada_c);  
   if  (2*L>=CE)
     lada_min= asin ((r-sqrt(r^2-[(1+lambda^2)*(r^2-lambda^2*R^2)]))/(R*(1+lambda^2)));   
   else
     lada_min= asin ((sqrt(4*L^2+r^2-(lambda*R)^2)*r-[2*L*lambda*R])/(4*L^2+r^2));    
   end
   

   counter_lada=1;  
 
    for lada= linspace(lada_min, lada_max,points) 
     
        nmin= ([r*cos(lada)]- [lambda*R])/(sin (lada)); 
        counter_n=1;
            for n=linspace(nmin, 2*L,points)  
           
                DF= sqrt(((1+lambda^2)*R^2)-r^2); 
                
                if (n>= DF)
                    zeta= atan ((sqrt((R*sin(lada))^2-(r-lambda*R*cos(lada))^2))/([r*cos(lada)]-[lambda*R]));
                else
                    zeta= acos (([r*cos(lada)]-[lambda*R])/(n*sin(lada)));  
                end

                f_3(counter_n)= zeta/(pi*L);
                if (isnan (zeta)==1)
                  f_3(counter_n)= 0;  
                end
                
                counter_n=counter_n+1;   

            end

            int3= ((2*L-nmin)/(2*m))*(f_3(1)+sum(2*f_3(2:end-1))+f_3(end)); 
            f_2(counter_lada)=sin(lada)*int3;
            counter_lada=counter_lada+1;
    
    end

    int2= ((lada_max-lada_min)/(2*m))* (f_2(1)+sum(2*f_2(2:end-1))+f_2(end));
    f_1(counter_r)=2*r*int2;
    counter_r=counter_r+1;        
end
N3b= ((R*sqrt(1+lambda^2) - lambda*R)/(2*m))* (f_1(1)+sum(2*f_1(2:end-1))+f_1(end));
N3b= real (N3b*lambda_t);

%%
fprintf('N1: %f\n', N1);
fprintf('N2: %f\n', N2);

fprintf('N3a: %f\n', N3a);
fprintf('N3b: %f\n', N3b);


fprintf('N4a: %f\n', N4a);
fprintf('N4b: %f\n', N4b);

fprintf('Sum all N : %f\n', N1+N2+N3a+N3b+N4a+N4b);



N_total=N1+N2+N4a+N4b+N3a+N3b;
disp ("Diffusivity ratio is equal to: ")
kappa=exp(-N_total)
disp ("****") 


