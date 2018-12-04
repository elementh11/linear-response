function [datab, datas, datads, datav] = getHamiltonianForG(file)
% -----------------------------------------------------------------------
% construct k and E resolved H for getting the Green's functions
% construct velocity operators
% results are saved for later integration
% -----------------------------------------------------------------------

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

%%%%%%%%%%%%%% preliminary stuff
vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;


load(hfile);
matrices = datah.matrices;
nrpts = datah.nrpts;
dim = datah.num_wann;

load(kfile);
kpoints = datak.kpoints;
nk = size(kpoints,1);



for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;
    
    %%%%%%%% construct the bulk HW90 and principal layers H00 and H01
    %%%%%%%% construct the velocity operators
    HW90 = zeros(dim);
    H00 = zeros(dim);
    H01 = zeros(dim);
    vx = zeros(dim);
    vy = zeros(dim);
    vz = zeros(dim);
    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
        ham = matrix.ham;
        hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
        vxcontr = 1i * realdisp(1) * hcontr;
        vycontr = 1i * realdisp(2) * hcontr;
        vzcontr = 1i * realdisp(3) * hcontr;
        HW90 = HW90 + hcontr;
        vx = vx + vxcontr;
        vy = vy + vycontr;
        vz = vz + vzcontr;
        if delta(3) == 0
            H00 = H00 + hcontr;
        end
        if delta(3) == 1
            H01 = H01 + hcontr;
        end
    end
    
    %%%%%%% get H iteratively
    ec = 0;
    for energy = energylist
        ec = ec + 1;
        eps = energy * eye(dim);
        [ H ] = iterativeGreensFunction(H00, H01, eps);
        
        ham_bulk{kc,ec} = H.b;
        ham_surf{kc,ec} = H.s;
        ham_dualsurf{kc,ec} = H.ds;
    end
    
    vop_x{kc} = vx;
    vop_y{kc} = vy;
    vop_z{kc} = vz;
    
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


