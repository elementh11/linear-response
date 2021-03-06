function data = opticalConductivity(file)
% ---- optical conductivity tensor with Kubo formula and iterative Green's function method

warning('off','all');
%%%%%%%%%%%%%% reading the input parameters
in = fopen(file,'r');
while(~feof(in))
    line = fgetl(in);
    if(contains(line,'=') && ~contains(line,'%'))
        variablestring = line(1:regexp(line,'=')-1);
        value = eval(line(regexp(line,'=')+1:end));
        variable = genvarname(line(1:regexp(line,'=')-1));
        eval([variable '= value;']);
        data.params.(variable) = value;
    end
end

%%%%%%%%% preliminary preparation
vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;

kb = 1.38e-23;  qelec = 1.6e-19; kt = kb / qelec * temp;
de = (elist(end) - elist(1)) / (length(elist) -1);

load(hfile);
matrices = datah.matrices;
nrpts = datah.nrpts;
dim = datah.num_wann;

load(kfile);
kpoints = datak.kpoints;
nk = size(kpoints,1);


% s_xx = zeros(nk, length(omegalist));
s_xy = zeros(nk, length(omegalist));
% s_xz = zeros(nk, length(omegalist));
s_yx = zeros(nk, length(omegalist));
% s_yy = zeros(nk, length(omegalist));
% s_yz = zeros(nk, length(omegalist));
% s_zx = zeros(nk, length(omegalist));
% s_zy = zeros(nk, length(omegalist));
% s_zz = zeros(nk, length(omegalist));

% sx = zeros(dim);
% sy = zeros(dim);
sz = zeros(dim);
for c = 1:dim/2
    upind = 2*c - 1;
    dnind = 2*c;
    sz(upind,upind) = +1;
    sz(dnind,dnind) = -1;
%     sx(upind,dnind) = 1;
%     sx(dnind,upind) = 1;
%     sy(upind,dnind) = -1i;
%     sy(dnind,upind) = +1i;
end



for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;
    
    %%%%%%%% construct the bulk HW90 and principal layers H00 and H01
    %%%%%%%% and the velocity operators
    HW90 = zeros(dim);
    H00 = zeros(dim);
    H01 = zeros(dim);
    vx = zeros(dim);
    vy = zeros(dim);
    jxz = zeros(dim);
    jyz = zeros(dim);
%     vz = zeros(dim);
    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
        ham = matrix.ham;
        hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
        vxcontr = 1i * realdisp(1) * hcontr;
        vycontr = 1i * realdisp(2) * hcontr;
%         vzcontr = 1i * realdisp(3) * hcontr;
        HW90 = HW90 + hcontr;
        vx = vx + vxcontr;
        vy = vy + vycontr;
%         vz = vz + vzcontr;
        if delta(3) == 0
            H00 = H00 + hcontr;
        end
        if delta(3) == 1
            H01 = H01 + hcontr;
        end
    end
    jxz = (vx * sz + sz * vx) / 2;
    jyz = (vy * sz + sz * vy) / 2;
    
    wc = 0;
    for omega = omegalist
        wc = wc + 1;
%         xx = 0;
        xy = 0;
%         xz = 0;
        yx = 0;
%         yy = 0;
%         yz = 0;
%         zx = 0;
%         zy = 0;
%         zz = 0;
%         ec = 0;
        for energy = elist
            ec = ec + 1;
            eps = energy * eye(dim);
            [H, HS, HSD] = getIterativeGreenFunction(H00, H01, eps);
            
            Gret = inv(eps - H + 1i*eta);
            GSret = inv(eps - HS + 1i*eta);
            GSDret = inv(eps - HSD + 1i*eta);
            Gadv = inv(eps - H - 1i*eta);
            GSadv = inv(eps - HS - 1i*eta);
            GSDadv = inv(eps - HSD - 1i*eta);
            fk = 1 / (exp((energy - ef)/kt) + 1);
            
            eps = (energy + omega) * eye(dim);
            [H, HS, HSD] = getIterativeGreenFunction(H00, H01, eps);
            
            Gret2 = inv(eps - H + 1i*eta);
            GSret2 = inv(eps - HS + 1i*eta);
            GSDret2 = inv(eps - HSD + 1i*eta);
            Gadv2 = inv(eps - H - 1i*eta);
            GSadv2 = inv(eps - HS - 1i*eta);
            GSDadv2 = inv(eps - HSD - 1i*eta);
            fk2 = 1 / (exp((energy + omega - ef)/kt) + 1);
            
%             contr_xx = fk * vx * Gret2 * vx * (Gret - Gadv) ...
%                 + fk2 * vx * (Gret2 - Gadv2) * vx * Gadv;
%             contr_xx(isnan(contr_xx)) = 0;
%             xx = xx + trace(contr_xx) * de;
            
            contr_xy = fk * vx * Gret2 * vy * (Gret - Gadv) ...
                + fk2 * vx * (Gret2 - Gadv2) * vy * Gadv;
            contr_xy(isnan(contr_xy)) = 0;
            xy = xy + trace(contr_xy) * de;
            
%             contr_xz = fk * vx * Gret2 * vz * (Gret - Gadv) ...
%                 + fk2 * vx * (Gret2 - Gadv2) * vz * Gadv;
%             contr_xz(isnan(contr_xz)) = 0;
%             xz = xz + trace(contr_xz) * de;
            %--------------------------------
            contr_yx = fk * vy * Gret2 * vx * (Gret - Gadv) ...
                + fk2 * vy * (Gret2 - Gadv2) * vx * Gadv;
            contr_yx(isnan(contr_yx)) = 0;
            yx = yx + trace(contr_yx) * de;
            
%             contr_yy = fk * vy * Gret2 * vy * (Gret - Gadv) ...
%                 + fk2 * vy * (Gret2 - Gadv2) * vy * Gadv;
%             contr_yy(isnan(contr_yy)) = 0;
%             yy = yy + trace(contr_yy) * de;
%             
%             contr_yz = fk * vy * Gret2 * vz * (Gret - Gadv) ...
%                 + fk2 * vy * (Gret2 - Gadv2) * vz * Gadv;
%             contr_yz(isnan(contr_yz)) = 0;
%             yz = yz + trace(contr_yz) * de;
%             %----------------------------------
%             contr_zx = fk * vz * Gret2 * vx * (Gret - Gadv) ...
%                 + fk2 * vz * (Gret2 - Gadv2) * vx * Gadv;
%             contr_zx(isnan(contr_zx)) = 0;
%             zx = zx + trace(contr_zx) * de;
%             
%             contr_zy = fk * vz * Gret2 * vy * (Gret - Gadv) ...
%                 + fk2 * vz * (Gret2 - Gadv2) * vy * Gadv;
%             contr_zy(isnan(contr_zy)) = 0;
%             zy = zy + trace(contr_zy) * de;
%             
%             contr_zz = fk * vz * Gret2 * vz * (Gret - Gadv) ...
%                 + fk2 * vz * (Gret2 - Gadv2) * vz * Gadv;
%             contr_zz(isnan(contr_zz)) = 0;
%             zz = zz + trace(contr_zz) * de;
        end
%         s_xx(kc, wc) = xx;
        s_xy(kc, wc) = xy;
%         s_xz(kc, wc) = xz;
        s_yx(kc, wc) = yx;
%         s_yy(kc, wc) = yy;
%         s_yz(kc, wc) = yz;
%         s_zx(kc, wc) = zx;
%         s_zy(kc, wc) = zy;
%         s_zz(kc, wc) = zz;
    end
    
    %%%%% insert progress indicator
    if mod(kc, 50)==0
        prog = fopen('progress.txt','a');
        fprintf(prog,'calculating for k=%d/%d \n',kc,nk);
        fclose(prog);
    end
end

% data.cond_xx = s_xx;
data.cond_xy = s_xy;
% data.cond_xz = s_xz;
data.cond_yx = s_yx;
% data.cond_yy = s_yy;
% data.cond_yz = s_yz;
% data.cond_zx = s_zx;
% data.cond_zy = s_zy;
% data.cond_zz = s_zz;
save(outputfile,'data','-v7.3');


