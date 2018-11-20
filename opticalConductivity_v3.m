function [data1, data2] = opticalConductivity_v3(file)
% --- !!! restricted to -xy and -yx components for now !!!
% v3 feats: ---> for each energy, get G and store it ONCE,
%                then for each \omega just read G to do energy integration

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
        %data1.params.(variable) = value;
    end
end

%%%%%%%%% preliminary stuff
vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;

kb = 1.38e-23;  qelec = 1.6e-19; kt = kb / qelec * abs_temp;
ecmax = length(energylist);
de = (energylist(end) - energylist(1)) / (ecmax -1);
step = w_to_e_ratio * de;
nw = length(0: step : omegamax);

load(hfile);
matrices = datah.matrices;
nrpts = datah.nrpts;
dim = datah.num_wann;

load(kfile);
kpoints = datak.kpoints;
nk = size(kpoints,1);

for region = {'b','s','ds'}
    s_xy.(region{1}) = zeros(nk, nw);
    s_yx.(region{1}) = zeros(nk, nw);
    spin_xy.(region{1}) = zeros(nk, nw);
    spin_yx.(region{1}) = zeros(nk, nw);
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
    
    for i = 1 : ecmax
        for region = {'b','s','ds'}
            GR(i).(region{1}) = zeros(dim);
            GA(i).(region{1}) = zeros(dim);
        end
    end
    
    ec = 0;
    for energy = energylist
        ec = ec + 1;
        eps = energy * eye(dim);
        [ H ] = getIterativeGreenFunction(H00, H01, eps);
        
        for region = {'b','s','ds'}
            greenret.(region{1}) = inv(eps - H.(region{1}) + 1i*eta);
            greenret.(region{1})(isnan(greenret.(region{1}))) = 0;
            greenadv.(region{1}) = inv(eps - H.(region{1}) - 1i*eta);
            greenadv.(region{1})(isnan(greenadv.(region{1}))) = 0;
        end
        GR(ec) = greenret;
        GA(ec) = greenadv;
    end
    
    wc = 0;
    for omega = 0: step : omegamax
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
            if ec + wc - 1 <= ecmax
                fk1 = 1 / (exp((energy - ef)/kt) + 1);
                fk2 = 1 / (exp((energy + omega - ef)/kt) + 1);
                for region = {'b','s','ds'}
                    xy_contr.(region{1}) = fk1 * vx * GR(ec+wc-1).(region{1}) * vy * (GR(ec).(region{1}) - GA(ec).(region{1})) ...
                        + fk2 * vx * (GR(ec+wc-1).(region{1}) - GA(ec+wc-1).(region{1})) * vy * GA(ec).(region{1});
                    xy_cond.(region{1}) = xy_cond.(region{1}) + trace(xy_contr.(region{1})) * de;
                    
                    yx_contr.(region{1}) = fk1 * vy * GR(ec+wc-1).(region{1}) * vx * (GR(ec).(region{1}) - GA(ec).(region{1})) ...
                        + fk2 * vy * (GR(ec+wc-1).(region{1}) - GA(ec+wc-1).(region{1})) * vx * GA(ec).(region{1});
                    yx_cond.(region{1}) = yx_cond.(region{1}) + trace(yx_contr.(region{1})) * de;
                    
                    spin_xy_contr.(region{1}) = fk1 * jxz * GR(ec+wc-1).(region{1}) * vy * (GR(ec).(region{1}) - GA(ec).(region{1})) ...
                        + fk2 * jxz * (GR(ec+wc-1).(region{1}) - GA(ec+wc-1).(region{1})) * vy * GA(ec).(region{1});
                    spin_xy_cond.(region{1}) = spin_xy_cond.(region{1}) + trace(spin_xy_contr.(region{1})) * de;
                    
                    spin_yx_contr.(region{1}) = fk1 * jyz * GR(ec+wc-1).(region{1}) * vx * (GR(ec).(region{1}) - GA(ec).(region{1})) ...
                        + fk2 * jyz * (GR(ec+wc-1).(region{1}) - GA(ec+wc-1).(region{1})) * vx * GA(ec).(region{1});
                    spin_yx_cond.(region{1}) = spin_yx_cond.(region{1}) + trace(spin_yx_contr.(region{1})) * de;
                end
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
    if mod(kc, prog_step)==0
        prog = fopen('progress.txt','a');
        fprintf(prog,'calculating for k=%d/%d \n',kc,nk);
        fclose(prog);
    end
end

for region = {'b','s','ds'}
    data1.sigma_xy.(region{1}) = s_xy.(region{1});
    data1.sigma_yx.(region{1}) = s_yx.(region{1});
    data2.spin_sigma_xy.(region{1}) = spin_xy.(region{1});
    data2.spin_sigma_yx.(region{1}) = spin_yx.(region{1});
end
save(outputfile1,'data1','-v7.3');
save(outputfile2,'data2','-v7.3');


