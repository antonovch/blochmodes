function [reff, phiT, teff] = calcEffectiveCoeffsSym(lambda,data,h,model,ld_cutoff)
% [reff, phiT, teff] = calcEffectiveCoeffsSym(kx,lambda,data,h,model,bw_cutoff)
% calculates effective coefficients, either 2 Bloch wave model with 2 BW as the dominant one,
% or 3 BW model with either BW 2 or 3 as the dominant in the case when
% Bloch wave reflection and transmission coefficient are different on the two sides
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
	tt1 =data.th1(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
	tt2 =data.th2(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));
    tt3 =data.th3(ld_cutoff(1):ld_cutoff(2)).'*ones(size(h));

	switch model
		case 4
			reff = calc_reff33(rr11,rr12,rr13,rr21,rr22,rr23,rr31,rr32,rr33,phi1,phi2,phi3);
			teff=[];
			phiT=real(phi3)+unwrap(angle(reff));
		case 4i
			reff = calc_reff32(rr11,rr12,rr13,rr21,rr22,rr23,rr31,rr32,rr33,phi1,phi2,phi3);
			teff=[];
			phiT=real(phi2)+unwrap(angle(reff));
		case 3
			[reff, teff] = calc_reff2(rr11,rr12,rr21,rr22,tt1,tt2,phi1,phi2);
			phiT=real(phi2)+unwrap(angle(reff));
		case 3i
			[reff, teff] = calc_reff2(rr11,rr13,rr31,rr33,tt1,tt3,phi1,phi3);
			phiT=real(phi3)+unwrap(angle(reff));
	end


end

function reff123 = calc_reff33(rr11,rr12,rr13,rr21,rr22,rr23,rr31,rr32,rr33,phi1,phi2,phi3)
R1=rr11.*exp(1i*phi1);R2=rr22.*exp(1i*phi2);R3=rr33.*exp(1i*phi3);
R21=rr21.*exp(1i*phi2);R23=rr23.*exp(1i*phi2);R12=rr12.*exp(1i*phi1);
R13=rr13.*exp(1i*phi1);R31=rr31.*exp(1i*phi3);R32=rr32.*exp(1i*phi3);
alpha1=1./(1-R1.^2);
R32_1=(R32+alpha1.*R1.*R12.*R31)./(1-alpha1.*R12.*R21);
R23_1=(R23+alpha1.*R1.*R13.*R21)./(1-alpha1.*R13.*R31);
delta12=alpha1.*R12.*R31./(1-alpha1.*R12.*R21);
delta13=alpha1.*R13.*R21./(1-alpha1.*R13.*R31);
Reff12=(R2+alpha1.*R1.*R12.*R21)./(1-alpha1.*R12.*R21);
Reff13=(R3+alpha1.*R1.*R13.*R31)./(1-alpha1.*R13.*R31);
alphaeff12=1./(1-Reff12.^2);
a=alphaeff12.*Reff12.*R32_1+alphaeff12.*delta12;
b=alphaeff12.*Reff12.*delta12+alphaeff12.*R32_1;

reff123=exp(-1i*phi3).*(Reff13+R23_1.*a+delta13.*b)./(1-R23_1.*b-delta13.*a);
end

function reff123 = calc_reff32(rr11,rr12,rr13,rr21,rr22,rr23,rr31,rr32,rr33,phi1,phi2,phi3)
R1=rr11.*exp(1i*phi1);R2=rr22.*exp(1i*phi2);R3=rr33.*exp(1i*phi3);
R21=rr21.*exp(1i*phi2);R23=rr23.*exp(1i*phi2);R12=rr12.*exp(1i*phi1);
R13=rr13.*exp(1i*phi1);R31=rr31.*exp(1i*phi3);R32=rr32.*exp(1i*phi3);
alpha1=1./(1-R1.^2);
R32_1=(R32+alpha1.*R1.*R12.*R31)./(1-alpha1.*R12.*R21);
R23_1=(R23+alpha1.*R1.*R13.*R21)./(1-alpha1.*R13.*R31);
delta12=alpha1.*R12.*R31./(1-alpha1.*R12.*R21);
delta13=alpha1.*R13.*R21./(1-alpha1.*R13.*R31);
Reff12=(R2+alpha1.*R1.*R12.*R21)./(1-alpha1.*R12.*R21);
Reff13=(R3+alpha1.*R1.*R13.*R31)./(1-alpha1.*R13.*R31);
alphaeff13=1./(1-Reff13.^2);
a=alphaeff13.*Reff13.*R23_1+alphaeff13.*delta13;
b=alphaeff13.*Reff13.*delta13+alphaeff13.*R23_1;

reff123=exp(-1i*phi2).*(Reff12+R32_1.*a+delta12.*b)./(1-R32_1.*b-delta12.*a);
end

function [reff, teff] = calc_reff2(rr11,rr12,rr21,rr22,tt1,tt2,phi1,phi2)

alpha=1./(1-rr11.^2.*exp(2*1i*phi1));
reff=(rr22+alpha.*rr11.*rr21.*rr12.*exp(2*1i*phi1))./(1-alpha.*rr21.*rr12.*exp(1i*(phi1+phi2)));
teff=tt2+tt1.*alpha.*rr21.*exp(1i*phi1).*(reff.*exp(1i*phi2)+rr11.*exp(1i.*phi1));

end