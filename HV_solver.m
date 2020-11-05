function [Am,Bm,Cm,Dm,Al,Bl,Cl,Dl,length]=HV_solver(text_path,library_path)
     warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames'); %Just for getting rid of a warning. Does not effect the
                                                            %output of the function
    raw_data=fopen(text_path);
    count=1;
    Data=cell(1,17);
    while count~= 29
        satir=fgetl(raw_data);
        if count==2
            Data{1,1}=satir;
        elseif count==4
            Data{2,1}=satir;
        elseif count==6
            Data{3,1}=satir;
        elseif count==8
            Data{4,1}=satir;
        elseif count==10
            Data{5,1}=satir;
        elseif count==12
            Data{6,1}=satir;
        elseif count==13
            Data{7,1}=satir;
        elseif count==15
            Data{8,1}=satir;
        elseif count==16
            Data{9,1}=satir;
        elseif count==18
            Data{10,1}=satir;
        elseif count==19
            Data{11,1}=satir;
        elseif count==21
            Data{12,1}=satir;
        elseif count==22
            Data{13,1}=satir;
        elseif count==24
            Data{14,1}=satir;
        elseif count==25
            Data{15,1}=satir;
        elseif count==27
            Data{16,1}=satir;
        elseif count==28
            Data{17,1}=satir;
        end
        count=count+1;
    end

    Ncircuit=str2double(Data{1,1});
    Nbundle=str2double(Data{2,1});
    dbundle=str2double(Data{3,1});
    length=str2double(Data{4,1}); %left as km for future use
    line_type=convertCharsToStrings(Data{5,1});
    D1AB=sqrt((str2double(Data(6,1))-str2double(Data(8,1)))^2 + (str2double(Data(7,1))-str2double(Data(9,1)))^2);
    D1AC=sqrt((str2double(Data(6,1))-str2double(Data(10,1)))^2 + (str2double(Data(7,1))-str2double(Data(11,1)))^2);
    D1BC=sqrt((str2double(Data(8,1))-str2double(Data(10,1)))^2 + (str2double(Data(9,1))-str2double(Data(11,1)))^2);

    if Ncircuit==1
        D2AB=-1;
        D2AC=-1;
        D2BC=-1;
    else
        D2AB=sqrt((str2double(Data(12,1))-str2double(Data(14,1)))^2 + (str2double(Data(13,1))-str2double(Data(15,1)))^2);
        D2AC=sqrt((str2double(Data(12,1))-str2double(Data(16,1)))^2 + (str2double(Data(13,1))-str2double(Data(17,1)))^2);
        D2BC=sqrt((str2double(Data(14,1))-str2double(Data(16,1)))^2 + (str2double(Data(15,1))-str2double(Data(17,1)))^2);
    end
    
    %This part below just checks the validity of line type
    library_table=readtable(library_path);
    istrue=0;
    while istrue==0
    for i = 1:15
        if strcmp(line_type,library_table.CodeWord(i,1))
            istrue=1;
        end
    end

    if istrue==0
        line_type=input("Undefined line type. Please try another one: ",'s');   %Please dont use "" while giving the new input.
    end
    end
   
       
    
    
    %phase 2 starts
    % First, find out where we are in the library.
    for i = 1:15
        if strcmp(line_type,library_table.CodeWord(i,1))
            rownum=i;
        end
    end
    
    
    % RESISTANCE
    Rac_per_km = table2array(library_table(rownum,7))/(1.60934);
    R = Rac_per_km/(Nbundle*Ncircuit);  
    
    
    
    %Calculate Ds_bundle for Nbundle=1 to 8. Convert feet to meter
    if Nbundle==1
        Ds_bundle=table2array(library_table(rownum,8))*0.3048;
        Ds_bundle_cap=table2array(library_table(rownum,8))*0.0254/2;
    
    elseif Nbundle==2
            Ds_bundle=sqrt(table2array(library_table(rownum,8))*0.3048*dbundle);
            Ds_bundle_cap=sqrt(table2array(library_table(rownum,5))*(0.0254/2)*dbundle);
      
    elseif Nbundle==3
            Ds_bundle=nthroot(table2array(library_table(rownum,8))*0.3048*(dbundle^2),3);
            Ds_bundle_cap=nthroot(table2array(library_table(rownum,5))*(0.0254/2)*(dbundle^2),3);
            
    elseif Nbundle==4
            Ds_bundle=nthroot(table2array(library_table(rownum,8))*0.3048*(dbundle^2)*(dbundle*sqrt(2)),4);
            Ds_bundle_cap=nthroot(table2array(library_table(rownum,5))*(0.0254/2)*(dbundle^2)*(dbundle*sqrt(2)),4);
            
    elseif Nbundle==5
            Ds_bundle=nthroot(table2array(library_table(rownum,8))*0.3048*(dbundle^2)*((dbundle*(1+sqrt(5))/2)^2),5); 
            Ds_bundle_cap=nthroot(table2array(library_table(rownum,5))*(0254/2)*(dbundle^2)*((dbundle*(1+sqrt(5))/2)^2),5);
            
   elseif Nbundle==6
            Ds_bundle=nthroot(table2array(library_table(rownum,8))*0.3048*(dbundle^2)*(2*dbundle)*((dbundle*sqrt(3))^2),6); 
            Ds_bundle_cap=nthroot(table2array(library_table(rownum,5))*(0.0254/2)*(dbundle^2)*(2*dbundle)*((dbundle*sqrt(3))^2),6);
            
   elseif Nbundle==7
       d=dbundle/(2*sin(pi/2/7));
       e=2*dbundle*cos(pi/7);
            Ds_bundle=nthroot(table2array(library_table(rownum,8))*0.3048*(dbundle^2)*(d^2)*(e^2),7); 
            Ds_bundle_cap=nthroot(table2array(library_table(rownum,5))*(0.0254/2)*(dbundle^2)*(d^2)*(e^2),7);  
            
   elseif Nbundle==8
       a=dbundle*sqrt(2);
       b=dbundle(1+2/sqrt(2));
       c=2*dbundle;
       Ds_bundle=nthroot(table2array(library_table(rownum,8))*0.3048*(dbundle^2)*(b^2)*(a^2)*c,8);
       Ds_bundle_cap=nthroot(table2array(library_table(rownum,5))*(0.0254/2)*(dbundle^2)*(b^2)*(a^2)*c,8);
    end
  
    
    
    
    %GMR and GMD values
    
    if Ncircuit==2
        
        %GMR in three different positions;
        dist_a1_a2=sqrt((str2double(Data(6,1))-str2double(Data(12,1)))^2 + (str2double(Data(7,1))-str2double(Data(13,1)))^2);
        GMR_a1_a2=sqrt(Ds_bundle*dist_a1_a2);
        GMR_a1_a2_cap=sqrt(Ds_bundle_cap*dist_a1_a2);
        
        dist_b1_b2=sqrt((str2double(Data(8,1))-str2double(Data(14,1)))^2 + (str2double(Data(9,1))-str2double(Data(15,1)))^2);
        GMR_b1_b2=sqrt(Ds_bundle*dist_b1_b2);
        GMR_b1_b2_cap=sqrt(Ds_bundle_cap*dist_b1_b2);
        
        dist_c1_c2=sqrt((str2double(Data(10,1))-str2double(Data(16,1)))^2 + (str2double(Data(11,1))-str2double(Data(17,1)))^2);
        GMR_c1_c2=sqrt(Ds_bundle*dist_c1_c2);
        GMR_c1_c2_cap=sqrt(Ds_bundle_cap*dist_c1_c2);
        
        GMR_p=nthroot((GMR_a1_a2*GMR_b1_b2*GMR_c1_c2),3);
        GMR_p_cap=nthroot((GMR_a1_a2_cap*GMR_b1_b2_cap*GMR_c1_c2_cap),3);
        
     
        %GMD timee!!
        % GMD-AB
        Da1_b2=sqrt((str2double(Data(6,1))-str2double(Data(14,1)))^2 + (str2double(Data(7,1))-str2double(Data(15,1)))^2);
        Da2_b1=sqrt((str2double(Data(12,1))-str2double(Data(8,1)))^2 + (str2double(Data(13,1))-str2double(Data(9,1)))^2);
        GMD_ab=nthroot(D1AB*Da1_b2*Da2_b1*D2AB,4);  
        % GMD-BC
        Db1_c2=sqrt((str2double(Data(8,1))-str2double(Data(16,1)))^2 + (str2double(Data(9,1))-str2double(Data(17,1)))^2);
        Db2_c1=sqrt((str2double(Data(14,1))-str2double(Data(10,1)))^2 + (str2double(Data(15,1))-str2double(Data(11,1)))^2);
        GMD_bc=nthroot(D1BC*Db1_c2*Db2_c1*D2BC,4);
        % GMD-CA
        Dc1_a2=sqrt((str2double(Data(10,1))-str2double(Data(12,1)))^2 + (str2double(Data(11,1))-str2double(Data(13,1)))^2);
        Dc2_a1=sqrt((str2double(Data(16,1))-str2double(Data(6,1)))^2 + (str2double(Data(17,1))-str2double(Data(7,1)))^2);
        GMD_cap=nthroot(D1AC*Dc1_a2*Dc2_a1*D2AC,4);
       
        GMD_2cct=nthroot((GMD_ab*GMD_bc*GMD_cap),3);
        
        %EFFECT OF EARTH FOR 2CCT
        ha1=2*(str2double(Data(7,1)));
        hb1=2*(str2double(Data(9,1)));
        hc1=2*(str2double(Data(11,1)));
        ha2=2*str2double(Data(13,1));
        hb2=2*str2double(Data(15,1));
        hc2=2*str2double(Data(17,1));
        ha1a2=sqrt((str2double(Data(6,1))-str2double(Data(12,1)))^2 + (str2double(Data(7,1))+str2double(Data(13,1)))^2);
        hb1b2=sqrt((str2double(Data(8,1))-str2double(Data(14,1)))^2 + (str2double(Data(9,1))+str2double(Data(15,1)))^2);
        hc1c2=sqrt((str2double(Data(10,1))-str2double(Data(16,1)))^2 + (str2double(Data(11,1))+str2double(Data(17,1)))^2);
        hden=nthroot(ha1*hb1*hc1*ha2*hb2*hc2*(ha1a2^2)*(hb1b2^2)*(hc1c2^2),12);
        
        ha1b1=sqrt((str2double(Data(6,1))-str2double(Data(8,1)))^2 + (str2double(Data(7,1))+str2double(Data(9,1)))^2);
        ha1c1=sqrt((str2double(Data(6,1))-str2double(Data(10,1)))^2 + (str2double(Data(7,1))+str2double(Data(11,1)))^2);
        hb1c1=sqrt((str2double(Data(8,1))-str2double(Data(10,1)))^2 + (str2double(Data(9,1))+str2double(Data(11,1)))^2);
        ha1b2=sqrt((str2double(Data(6,1))-str2double(Data(14,1)))^2 + (str2double(Data(7,1))+str2double(Data(15,1)))^2);
        ha1c2=sqrt((str2double(Data(16,1))-str2double(Data(6,1)))^2 + (str2double(Data(17,1))+str2double(Data(7,1)))^2);
        hb1c2=sqrt((str2double(Data(8,1))-str2double(Data(16,1)))^2 + (str2double(Data(9,1))+str2double(Data(17,1)))^2);
        hb1a2=sqrt((str2double(Data(12,1))-str2double(Data(8,1)))^2 + (str2double(Data(13,1))+str2double(Data(9,1)))^2);
        hc1a2=sqrt((str2double(Data(10,1))-str2double(Data(12,1)))^2 + (str2double(Data(11,1))+str2double(Data(13,1)))^2);
        hc1b2=sqrt((str2double(Data(14,1))-str2double(Data(10,1)))^2 + (str2double(Data(15,1))+str2double(Data(11,1)))^2);
        ha2b2=sqrt((str2double(Data(12,1))-str2double(Data(14,1)))^2 + (str2double(Data(13,1))+str2double(Data(15,1)))^2);
        ha2c2=sqrt((str2double(Data(12,1))-str2double(Data(16,1)))^2 + (str2double(Data(13,1))+str2double(Data(17,1)))^2);
        hb2c2=sqrt((str2double(Data(14,1))-str2double(Data(16,1)))^2 + (str2double(Data(15,1))+str2double(Data(17,1)))^2);
        hnum=nthroot(ha1b1*ha1c1*hb1c1*ha1b2*ha1c2*hb1c2*hb1a2*hc1a2*hc1b2*ha2b2*ha2c2*hb2c2,12);
       
        L=2*(10^(-7))*log(GMD_2cct/GMR_p);
        C=2*pi*8.8541878128*(10^-12)/(log(GMD_2cct/GMR_p_cap)-log(hnum/hden));
        
        
    elseif Ncircuit==1
        
        GMD_1cct=nthroot((D1AB*D1AC*D1BC),3);
        
        %EFFECT OF EARTH FOR 1CCT
        ha=2*(str2double(Data(7,1)));
        hb=2*(str2double(Data(9,1)));
        hc=2*(str2double(Data(11,1)));
        hden=nthroot(ha*hb*hc,3);
        
        hab=sqrt((str2double(Data(6,1))-str2double(Data(8,1)))^2 + (str2double(Data(7,1))+str2double(Data(9,1)))^2);
        hac=sqrt((str2double(Data(6,1))-str2double(Data(10,1)))^2 + (str2double(Data(7,1))+str2double(Data(11,1)))^2);
        hbc=sqrt((str2double(Data(8,1))-str2double(Data(10,1)))^2 + (str2double(Data(9,1))+str2double(Data(11,1)))^2);
        hnum=nthroot(hab*hac*hbc,3);
        
        C=2*pi*8.85*(10^-12)/(log(GMD_1cct/Ds_bundle_cap)-log(hnum/hden));
        L=2*(10^(-7))*log(GMD_1cct/Ds_bundle);
    end
    
    
    X=2*pi*50*L*1000;
    B=2*pi*50*C*1000;
    XBonus=-1;
    BBonus=-1;
    
        %phase 3 starts here
    
    series_impedance=(R+1i*X)*length; %ohm/km * km
    shunt_admittance = 1i*B*length; % mho/km * km
    
    %MEDIUM LINE MODEL
    Am = (1+(series_impedance*shunt_admittance)/2);
    Dm=Am; %unitless
    Bm = series_impedance; %ohm
    Cm = shunt_admittance*(1+(series_impedance*shunt_admittance)/4); %mho
    
    %LONG LINE MODEL
    Zo=sqrt(series_impedance/shunt_admittance); %characteristic impedance
    gamma=sqrt(series_impedance*shunt_admittance); %propogation const.
    
    Al=cosh(gamma);
    Bl=Zo*sinh(gamma);
    Cl=sinh(gamma)/Zo;
    Dl=cosh(gamma);
end

