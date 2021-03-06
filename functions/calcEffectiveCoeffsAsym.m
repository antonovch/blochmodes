function [reff, reff2, phiT, teff, teff2] = calcEffectiveCoeffsAsym(lambda,data,data2,h,model,ld_cutoff)
% [reff, reff2, phiT, teff, teff2] = calcEffectiveCoeffsAsym(kx,lambda,data,data2,h,model,bw_cutoff)
% calculates effective coefficients, either 2 Bloch wave model with 2 BW as the dominant one,
% or 3 BW model with either BW 2 or 3 as the dominant
% ld_cutoff vector of length 2 with indices of start/stop for the wavelength
% Not implemented: effective transmission for 3 BW model

	phi1=real((2*pi./lambda(ld_cutoff(1):ld_cutoff(2)).*data.n_eff(1,ld_cutoff(1):ld_cutoff(2))).'*h);
	phi2=real((2*pi./lambda(ld_cutoff(1):ld_cutoff(2)).*data.n_eff(2,ld_cutoff(1):ld_cutoff(2))).'*h);
	phi3=real((2*pi./lambda(ld_cutoff(1):ld_cutoff(2)).*data.n_eff(3,ld_cutoff(1):ld_cutoff(2))).'*h);
	rr11=data.r11(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr22=data.r22(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr12=data.r12(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr13=data.r13(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr23=data.r23(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr33=data.r33(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr31=data.r31(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr32=data.r32(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr21=data.r21(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	tt1=data.th1(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	tt2=data.th2(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr11_2=data2.r11(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr22_2=data2.r22(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr12_2=data2.r12(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr13_2=data2.r13(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr23_2=data2.r23(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr33_2=data2.r33(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr31_2=data2.r31(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr32_2=data2.r32(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	rr21_2=data2.r21(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	tt1_2=data2.th1(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	tt2_2=data2.th2(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));

	switch model
		case 4
			[reff, reff2, teff, teff2] = calc_reff33(rr11,rr11_2,rr12,rr12_2,...
            rr13,rr13_2,rr21,rr21_2,rr22,rr22_2,rr23,rr23_2,rr31,rr31_2,...
            rr32,rr32_2,rr33,rr33_2,phi1,phi2,phi3);
			phiT=real(phi3)+(unwrap(angle(reff))+unwrap(angle(reff2)))/2;
		case 4i
			[reff, reff2, teff, teff2] = calc_reff32(rr11,rr11_2,rr12,rr12_2,...
            rr13,rr13_2,rr21,rr21_2,rr22,rr22_2,rr23,rr23_2,rr31,rr31_2,...
            rr32,rr32_2,rr33,rr33_2,phi1,phi2,phi3);
			phiT=real(phi2)+(unwrap(angle(reff))+unwrap(angle(reff2)))/2;
		case 3
			[reff, reff2, teff, teff2] = calc_reff2(rr11,rr11_2,rr12,rr12_2,...
                rr21,rr21_2,rr22,rr22_2,tt1,tt1_2,tt2,tt2_2,phi1,phi2);
			phiT=real(phi2)+(unwrap(angle(reff))+unwrap(angle(reff2)))/2;
		case 3i
			[reff, reff2, teff, teff2] = calc_reff2(rr11,rr11_2,rr13,rr13_2,...
                rr31,rr31_2,rr33,rr33_2,tt1,tt1_2,tt3,tt3_2,phi1,phi3);
			phiT=real(phi3)+(unwrap(angle(reff))+unwrap(angle(reff2)))/2;
	end


end

function [reff123, reff123_2, teff, teff2] = calc_reff33(rr11,rr11_2,rr12,rr12_2,...
            rr13,rr13_2,rr21,rr21_2,rr22,rr22_2,rr23,rr23_2,rr31,rr31_2,...
            rr32,rr32_2,rr33,rr33_2,phi1,phi2,phi3)
        
R1=rr11.*exp(1i*phi1);R2=rr22.*exp(1i*phi2);R3=rr33.*exp(1i*phi3);
R21=rr21.*exp(1i*phi2);R23=rr23.*exp(1i*phi2);R12=rr12.*exp(1i*phi1);
R13=rr13.*exp(1i*phi1);R31=rr31.*exp(1i*phi3);R32=rr32.*exp(1i*phi3);
R1_2=rr11_2.*exp(1i*phi1);R2_2=rr22_2.*exp(1i*phi2);R3_2=rr33_2.*exp(1i*phi3);
R21_2=rr21_2.*exp(1i*phi2);R23_2=rr23_2.*exp(1i*phi2);R12_2=rr12_2.*exp(1i*phi1);
R13_2=rr13_2.*exp(1i*phi1);R31_2=rr31_2.*exp(1i*phi3);R32_2=rr32_2.*exp(1i*phi3);
alpha1=1./(1-R1_2.*R1);

R32_1=(R32+alpha1.*R1_2.*R12.*R31)./(1-alpha1.*R12.*R21_2);
R32_1_2=(R32_2+alpha1.*R1.*R12_2.*R31_2)./(1-alpha1.*R12_2.*R21);
R23_1=(R23+alpha1.*R1_2.*R13.*R21)./(1-alpha1.*R13.*R31_2);
R23_1_2=(R23_2+alpha1.*R1.*R13_2.*R21_2)./(1-alpha1.*R13_2.*R31);
delta12=alpha1.*R12.*R31_2./(1-alpha1.*R12.*R21_2);
delta12_2=alpha1.*R12_2.*R31./(1-alpha1.*R12_2.*R21);
delta13=alpha1.*R13.*R21_2./(1-alpha1.*R13.*R31_2);
delta13_2=alpha1.*R13_2.*R21./(1-alpha1.*R13_2.*R31);
Reff12=(R2+alpha1.*R1_2.*R12.*R21)./(1-alpha1.*R12.*R21_2);
Reff12_2=(R2_2+alpha1.*R1.*R12_2.*R21_2)./(1-alpha1.*R12_2.*R21);
Reff13=(R3+alpha1.*R1_2.*R13.*R31)./(1-alpha1.*R13.*R31_2);
Reff13_2=(R3_2+alpha1.*R1.*R13_2.*R31_2)./(1-alpha1.*R13_2.*R31);
alphaeff12=1./(1-Reff12.*Reff12_2);
a=alphaeff12.*Reff12_2.*R32_1+alphaeff12.*delta12_2;
a_2=alphaeff12.*Reff12.*R32_1_2+alphaeff12.*delta12;
b=alphaeff12.*Reff12_2.*delta12+alphaeff12.*R32_1_2;
b_2=alphaeff12.*Reff12.*delta12_2+alphaeff12.*R32_1;

reff123=exp(-1i*phi3).*(Reff13+R23_1.*a+delta13.*b_2)./(1-R23_1.*b-delta13.*a_2);
reff123_2=exp(-1i*phi3).*(Reff13_2+R23_1_2.*a_2+delta13_2.*b)./(1-R23_1_2.*b_2-delta13_2.*a);
teff=[];teff2=[];
end

function [reff123, reff123_2, teff, teff2] = calc_reff32(rr11,rr11_2,rr12,rr12_2,...
            rr13,rr13_2,rr21,rr21_2,rr22,rr22_2,rr23,rr23_2,rr31,rr31_2,...
            rr32,rr32_2,rr33,rr33_2,phi1,phi2,phi3)
        
R1=rr11.*exp(1i*phi1);R2=rr22.*exp(1i*phi2);R3=rr33.*exp(1i*phi3);
R21=rr21.*exp(1i*phi2);R23=rr23.*exp(1i*phi2);R12=rr12.*exp(1i*phi1);
R13=rr13.*exp(1i*phi1);R31=rr31.*exp(1i*phi3);R32=rr32.*exp(1i*phi3);
R1_2=rr11_2.*exp(1i*phi1);R2_2=rr22_2.*exp(1i*phi2);R3_2=rr33_2.*exp(1i*phi3);
R21_2=rr21_2.*exp(1i*phi2);R23_2=rr23_2.*exp(1i*phi2);R12_2=rr12_2.*exp(1i*phi1);
R13_2=rr13_2.*exp(1i*phi1);R31_2=rr31_2.*exp(1i*phi3);R32_2=rr32_2.*exp(1i*phi3);
alpha1=1./(1-R1_2.*R1);

R32_1=(R32+alpha1.*R1_2.*R12.*R31)./(1-alpha1.*R12.*R21_2);
R32_1_2=(R32_2+alpha1.*R1.*R12_2.*R31_2)./(1-alpha1.*R12_2.*R21);
R23_1=(R23+alpha1.*R1_2.*R13.*R21)./(1-alpha1.*R13.*R31_2);
R23_1_2=(R23_2+alpha1.*R1.*R13_2.*R21_2)./(1-alpha1.*R13_2.*R31);
delta12=alpha1.*R12.*R31_2./(1-alpha1.*R12.*R21_2);
delta12_2=alpha1.*R12_2.*R31./(1-alpha1.*R12_2.*R21);
delta13=alpha1.*R13.*R21_2./(1-alpha1.*R13.*R31_2);
delta13_2=alpha1.*R13_2.*R21./(1-alpha1.*R13_2.*R31);
Reff12=(R2+alpha1.*R1_2.*R12.*R21)./(1-alpha1.*R12.*R21_2);
Reff12_2=(R2_2+alpha1.*R1.*R12_2.*R21_2)./(1-alpha1.*R12_2.*R21);
Reff13=(R3+alpha1.*R1_2.*R13.*R31)./(1-alpha1.*R13.*R31_2);
Reff13_2=(R3_2+alpha1.*R1.*R13_2.*R31_2)./(1-alpha1.*R13_2.*R31);
alphaeff13=1./(1-Reff13.*Reff13_2);
a=alphaeff13.*Reff13_2.*R23_1+alphaeff13.*delta13_2;
a_2=alphaeff13.*Reff13.*R23_1_2+alphaeff13.*delta13;
b=alphaeff13.*Reff13_2.*delta13+alphaeff13.*R23_1_2;
b_2=alphaeff13.*Reff13.*delta13_2+alphaeff13.*R23_1;

reff123=exp(-1i*phi2).*(Reff12+R32_1.*a+delta12.*b_2)./(1-R32_1.*b-delta12.*a_2);
reff123_2=exp(-1i*phi2).*(Reff12_2+R32_1_2.*a_2+delta12_2.*b)./(1-R32_1_2.*b_2-delta12_2.*a);
teff=[];teff2=[];
end

function [reff, reff2, teff, teff2] = calc_reff2(rr11,rr11_2,rr12,rr12_2,rr21,rr21_2,rr22,rr22_2,tt1,tt1_2,tt2,tt2_2,phi1,phi2)

alpha=1./(1-rr11.*rr11_2.*exp(2*1i*phi1)); % this mixing is why we cannot do two individual runs of the model
reff=(rr22+alpha.*rr11_2.*rr21.*rr12.*exp(2*1i*phi1))./(1-alpha.*rr21.*rr12_2.*exp(1i*(phi1+phi2)));
reff2=(rr22_2+alpha.*rr11.*rr21_2.*rr12_2.*exp(2*1i*phi1))./(1-alpha.*rr21_2.*rr12.*exp(1i*(phi1+phi2)));
teff=[];teff2=[];

end