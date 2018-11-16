function data1 = opticalConductivity_v2(file)
% ---- optical conductivity tensor with Kubo formula and iterative Green's function method
% --- !!! restricted to -xy and -yx components for now !!!

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
        data1.params.(variable) = value;
    end
end

%%%%%%%%% preliminary stuff
vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;

kb = 1.38e-23;  qelec = 1.6e-19; kt = kb / qelec * temp;
de = (energylist(end) - energylist(1)) / (length(energylist) -1);

load(hfile);
matrices = datah.matrices;
nrpts = datah.nrpts;
dim = datah.num_wann;

load(kfile);
kpoints = datak.kpoints;
nk = size(kpoints,1);

for region = {'b','s','ds'}
    s_xy.(region{1}) = zeros(nk, length(omegalist));
    s_yx.(region{1}) = zeros(nk, length(omegalist));
    spin_xy.(region{1}) = zeros(nk, length(omegalist));
    spin_yx.(region{1}) = zeros(nk, length(omegalist));
end

sz = zeros(dim);
for c = 1:dim/2
    upind = 2*c - 1;
    dnind = 2*c;
    sz(upind,upind) = +1;
    sz(dnind,dnind) = -1;
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
    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
        ham = matrix.ham;
        hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
        vxcontr = 1i * realdisp(1) * hcontr;
        vycontr = 1i * realdisp(2) * hcontr;
        HW90 = HW90 + hcontr;
        vx = vx + vxcontr;
        vy = vy + vycontr;
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
        for region = {'b','s','ds'}
            xy_cond.(region{1}) = 0;
            yx_cond.(region{1}) = 0;
            spin_xy_cond.(region{1}) = 0;
            spin_yx_cond.(region{1}) = 0;
        end
        ec = 0;
        for energy = energylist
            ec = ec + 1;
            fk1 = 1 / (exp((energy - ef)/kt) + 1);
            fk2 = 1 / (exp((energy + omega - ef)/kt) + 1);
            eps1 = energy * eye(dim);
            eps2 = (energy + omega) * eye(dim);
            [ H1 ] = getIterativeGreenFunction(H00, H01, eps1);
            [ H2 ] = getIterativeGreenFunction(H00, H01, eps2);
            
            for region = {'b','s','ds'}
                G1R.(region{1}) = inv(eps1 - H1.(region{1}) + 1i*eta);
                G1A.(region{1}) = inv(eps1 - H1.(region{1}) - 1i*eta);
                G2R.(region{1}) = inv(eps2 - H2.(region{1}) + 1i*eta);
                G2A.(region{1}) = inv(eps2 - H2.(region{1}) - 1i*eta);
                
                %------ Kubo formula in terms of Green's functions
                xy_contr.(region{1}) = fk1 * vx * G2R.(region{1}) * vy * (G1R.(region{1}) - G1A.(region{1})) ...
                    + fk2 * vx * (G2R.(region{1}) - G2A.(region{1})) * vy * G1A.(region{1});
                %------ correction for singularities
                xy_contr.(region{1})(isnan(xy_contr.(region{1}))) = 0;
                %------ perform energy integration
                xy_cond.(region{1}) = xy_cond.(region{1}) + trace(xy_contr.(region{1})) * de;
                
                
                yx_contr.(region{1}) = fk1 * vy * G2R.(region{1}) * vx * (G1R.(region{1}) - G1A.(region{1})) ...
                    + fk2 * vy * (G2R.(region{1}) - G2A.(region{1})) * vx * G1A.(region{1});
                yx_contr.(region{1})(isnan(yx_contr.(region{1}))) = 0;
                yx_cond.(region{1}) = yx_cond.(region{1}) + trace(yx_contr.(region{1})) * de;
                
                spin_xy_contr.(region{1}) = fk1 * jxz * G2R.(region{1}) * vy * (G1R.(region{1}) - G1A.(region{1})) ...
                    + fk2 * jxz * (G2R.(region{1}) - G2A.(region{1})) * vy * G1A.(region{1});
                spin_xy_contr.(region{1})(isnan(spin_xy_contr.(region{1}))) = 0;
                spin_xy_cond.(region{1}) = spin_xy_cond.(region{1}) + trace(spin_xy_contr.(region{1})) * de;
                
                spin_yx_contr.(region{1}) = fk1 * jyz * G2R.(region{1}) * vx * (G1R.(region{1}) - G1A.(region{1})) ...
                    + fk2 * jyz * (G2R.(region{1}) - G2A.(region{1})) * vx * G1A.(region{1});
                spin_yx_contr.(region{1})(isnan(spin_yx_contr.(region{1}))) = 0;
                spin_yx_cond.(region{1}) = spin_yx_cond.(region{1}) + trace(spin_yx_contr.(region{1})) * de;
            end
            
        end
        for region = {'b','s','ds'}
            s_xy.(region{1})(kc, wc) = xy_cond.(region{1});
            s_yx.(region{1})(kc, wc) = yx_cond.(region{1});
            spin_xy.(region{1})(kc, wc) = spin_xy_cond.(region{1});
            spin_yx.(region{1})(kc, wc) = spin_yx_cond.(region{1});
        end
    end
    
    %%%%%%%%%%%%%%%% progress indicator
    if mod(kc, 50)==0
        prog = fopen('progress.txt','a');
        fprintf(prog,'calculating for k=%d/%d \n',kc,nk);
        fclose(prog);
    end
end

for region = {'b','s','ds'}
    data1.sigma_xy.(region{1}) = s_xy.(region{1});
    data1.sigma_yx.(region{1}) = s_yx.(region{1});
    data1.spin_sigma_xy.(region{1}) = spin_xy.(region{1});
    data1.spin_sigma_yx.(region{1}) = spin_yx.(region{1});
end
save(outputfile1,'data1','-v7.3');


