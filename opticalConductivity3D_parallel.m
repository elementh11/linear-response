function [datab, datas, datads, datav] = opticalConductivity3D_parallel(file)
% -----------------------------------------------------------------------
% read H and v operators from "getHamiltonianForG_parallel.m" output
% get optical sigma tensor components
% save \omega and k resolved sigma
% -----------------------------------------------------------------------
%parpool(4);
warning('off','all');
%tic
%%%%%%%%%%%%%% reading the input parameters
in = fopen(file,'r');
while(~feof(in))
    line = fgetl(in);
    if(contains(line,'=') && ~contains(line,'%'))
        variablestring = line(1:regexp(line,'=')-1);
        value = eval(line(regexp(line,'=')+1:end));
        variable = matlab.lang.makeValidName(line(1:regexp(line,'=')-1));
        eval([variable '= value;']);
        %data1.params.(variable) = value;
    end
end

%%%%%%%%%%%%%% preliminary stuff
ef=ef;energylist=energylist;energyval=energyval;omegaval=omegaval;
kb = 1.38e-23;  qelec = 1.6e-19; kt = kb / qelec * abs_temp;

load(kfile);
kpoints = datak.kpoints;
nk = size(kpoints,1);
nw = length(omegaval);

for r = {'b','s','ds'}
    for cpt = {'xx','xy','xz','yx','yy','yz','zx','zy','zz'}
        sigma_kw.(r{1}).(cpt{1}) = zeros(nw, nk);
        sigma_w.(r{1}).(cpt{1}) = zeros(nw, 1);
    end
end

wc = 0;
for omega = omegaval
    wc = wc + 1;
    windex = wc - 1;
    for r = {'b','s','ds'}
        for cpt = {'xx','xy','xz','yx','yy','yz','zx','zy','zz'}
            sigma_k.(r{1}).(cpt{1}) = zeros(nk, 1);
        end
    end
    
    for kc = 1:nk
        warning('off','all');
        
        for r = {'b','s','ds'}
            for cpt = {'xx','xy','xz','yx','yy','yz','zx','zy','zz'}
                sigma_temp.(r{1}).(cpt{1}) = 0;
            end
        end
        
        ec = 0;
        for energy = energyval
            ec = ec + 1;
            eps = energy * eye(dim);
            fk1 = 1 / (exp((energy - ef)/kt) + 1);
            fk2 = 1 / (exp((energy + omega - ef)/kt) + 1);
            
        for r = {'b','s','ds'}
            
            GR1.(r{1}) = inv(eps -ham.(r{1}) + 1i*eta);
            
            
            
            
            
            for cpt = {'xx','xy','xz','yx','yy','yz','zx','zy','zz'}
                contr.(r{1}).(cpt{1}) = trace(v);
                sigma_temp.(r{1}).(cpt{1}) = sigma_temp.(r{1}).(cpt{1}) + contr.(r{1}).(cpt{1});
            end
        end
            
            
            
            
            
            
            
    
    
    
    
    
    
    %%%%%%% get H iteratively
    ec = 0;
    for energy = energylist
        ec = ec + 1;
        eps = energy * eye(dim);
        [ H ] = iterativeGreensFunction(H00, H01, eps);
        
        ham_b(ec, :, :) = H.b;
        ham_s(ec, :, :) = H.s;
        ham_ds(ec, :, :) = H.ds;
    end
    
    ham_bulk(kc, :, :, :) = ham_b;
    ham_surf(kc, :, :, :) = ham_s;
    ham_dualsurf(kc, :, :, :) = ham_ds;
    
    vop_x(kc, :, :) = vx;
    vop_y(kc, :, :) = vy;
    vop_z(kc, :, :) = vz;
    
    %%%%%%%%%%%%%%%% progress indicator
    fprintf('calculation for k = %d / %d done... \n',kc,nk);
    if mod(kc, prog_step)==0
        prog = fopen('progress.txt','a');
        fprintf(prog,'calculation for k = %d / %d done... \n',kc,nk);
        fclose(prog);
    end
end

datab.ham_bulk = ham_bulk;
datas.ham_surf = ham_surf;
datads.ham_dualsurf = ham_dualsurf;
datav.vop_x = vop_x;
datav.vop_y = vop_y;
datav.vop_z = vop_z;

save(fileb,'datab','-v7.3');
save(files,'datas','-v7.3');
save(fileds,'datads','-v7.3');
save(filev,'datav','-v7.3');
toc
delete(gcp('nocreate'));



