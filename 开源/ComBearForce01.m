function [F_xR,F_yR,resr_d] = ComBearForce01(w,r,R,N_b,t,cb,x,y,r_0,isfault)
wcage = w*(1-r/R)/2;
F_xR=0;
F_yR=0;
LD = 0.5e-3;
A = 2e-3;
rd = 9.5e-3;%滚动体
dout = 30e-3;
rd =rd/2;
deta = rd-sqrt(rd^2-(LD/2)^2);
qout =3*pi/2; 
resr_d = 0.0;

    for i=1:1:N_b
        r_d = 0;
        thea = wcage*t+2*pi/N_b*(i-1);
        qout = w*t;
        aaa = asin(sin(abs(qout-thea)));
        betat = asin(LD/(dout));
       
        if isfault == 1
            if abs(mod(abs(thea-qout),2*pi))<=betat
                r_d = deta;
            end
        else
            r_d = 0;
        end
        if r_d>0
           resr_d =  r_d;
        end
        
        inp = x*cos(thea)+y*sin(thea)-r_0-r_d;
        if (inp>0)
            c1 = cb*inp^1.5;
            F_xR=F_xR-c1*cos(thea);
            F_yR = F_yR-c1*sin(thea);
            
        else
            F_xR=F_xR+0;
            F_yR = F_yR+0;
            
        end
    end


end