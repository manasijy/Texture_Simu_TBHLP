function SlipSystem = SS_setFCC_fp_function()
CS = crystalSymmetry('-43m');
SlipSystem = struct('n',zeros(1,3),'b',zeros(1,3));

SlipSystem(1).n=[1,1,1]; SlipSystem(1).b =[0,-1,1];
SlipSystem(2).n=[1,1,1]; SlipSystem(2).b =[-1,0,1];
SlipSystem(3).n=[1,1,1]; SlipSystem(3).b =[-1,1,0];
SlipSystem(4).n=[-1,1,1]; SlipSystem(4).b =[1,1,0];
SlipSystem(5).n=[-1,1,1];SlipSystem(5).b =[1,0,1];
SlipSystem(6).n=[-1,1,1];SlipSystem(6).b =[0,-1,1];
SlipSystem(7).n=[1,-1,1]; SlipSystem(7).b =[1,1,0];
SlipSystem(8).n=[1,-1,1]; SlipSystem(8).b =[-1,0,1];
SlipSystem(9).n=[1,-1,1]; SlipSystem(9).b =[0,1,1];
SlipSystem(10).n=[1,1,-1];SlipSystem(10).b =[-1,1,0];
SlipSystem(11).n=[1,1,-1]; SlipSystem(11).b=[1,0,1];
SlipSystem(12).n=[1,1,-1]; SlipSystem(12).b =[0,1,1];
SlipSystem(13).n=[1,1,1]; SlipSystem(13).b =[-1,2,-1];
SlipSystem(14).n=[1,1,1]; SlipSystem(14).b =[-1,-1,2];
SlipSystem(15).n=[1,1,1]; SlipSystem(15).b =[2,-1,-1];
SlipSystem(16).n=[-1,1,1]; SlipSystem(16).b =[-1,1,-2];
SlipSystem(17).n=[-1,1,1];SlipSystem(17).b =[-1,-2,1];
SlipSystem(18).n=[-1,1,1];SlipSystem(18).b =[2,1,1];
SlipSystem(19).n=[1,-1,1]; SlipSystem(19).b =[1,-1,-2];
SlipSystem(20).n=[1,-1,1]; SlipSystem(20).b =[-2,-1,1];
SlipSystem(21).n=[1,-1,1]; SlipSystem(21).b =[1,2,1];
SlipSystem(22).n=[1,1,-1]; SlipSystem(22).b =[-2,1,-1];
SlipSystem(23).n=[1,1,-1];SlipSystem(23).b =[1,1,2];
SlipSystem(24).n=[1,1,-1];SlipSystem(24).b =[1,-2,-1];


SlipSystem(25).n=[1,1,1]; SlipSystem(25).b =[0,1,-1];
SlipSystem(26).n=[1,1,1]; SlipSystem(26).b =[1,0,-1];
SlipSystem(27).n=[1,1,1]; SlipSystem(27).b =[1,-1,0];
SlipSystem(28).n=[-1,1,1]; SlipSystem(28).b =[-1,-1,0];
SlipSystem(29).n=[-1,1,1];SlipSystem(29).b =[-1,0,-1];
SlipSystem(30).n=[-1,1,1];SlipSystem(30).b =[0,1,-1];
SlipSystem(31).n=[1,-1,1]; SlipSystem(31).b =[-1,-1,0];
SlipSystem(32).n=[1,-1,1]; SlipSystem(32).b =[1,0,-1];
SlipSystem(33).n=[1,-1,1]; SlipSystem(33).b =[0,-1,-1];
SlipSystem(34).n=[1,1,-1];SlipSystem(34).b =[1,-1,0];
SlipSystem(35).n=[1,1,-1]; SlipSystem(35).b=[-1,0,-1];
SlipSystem(36).n=[1,1,-1]; SlipSystem(36).b =[0,-1,-1];
SlipSystem(37).n=[1,1,1]; SlipSystem(37).b =[1,-2,1];
SlipSystem(38).n=[1,1,1]; SlipSystem(38).b =[1,1,-2];
SlipSystem(39).n=[1,1,1]; SlipSystem(39).b =[-2,1,1];
SlipSystem(40).n=[-1,1,1]; SlipSystem(40).b =[1,-1,2];
SlipSystem(41).n=[-1,1,1];SlipSystem(41).b =[1,2,-1];
SlipSystem(42).n=[-1,1,1];SlipSystem(42).b =[-2,-1,-1];
SlipSystem(43).n=[1,-1,1]; SlipSystem(43).b =[-1,1,2];
SlipSystem(44).n=[1,-1,1]; SlipSystem(44).b =[2,1,-1];
SlipSystem(45).n=[1,-1,1]; SlipSystem(45).b =[-1,-2,-1];
SlipSystem(46).n=[1,1,-1]; SlipSystem(46).b =[2,-1,1];
SlipSystem(47).n=[1,1,-1];SlipSystem(47).b =[-1,-1,-2];
SlipSystem(48).n=[1,1,-1];SlipSystem(48).b =[-1,2,1];

    for si = 1:1:48
            n = SlipSystem(si).n;
            b = SlipSystem(si).b;
            hkl = Miller(n(1),n(2),n(3),CS);
            uvw = Miller(b(1),b(2),b(3),CS);
            beta_ij = SchmidTensor(uvw,hkl,'generalized'); 
            m_ij = SchmidTensor(uvw,hkl,'symmetric');
            q_ij = SchmidTensor(uvw,hkl,'antisymmetric');
            SlipSystem(si).m = m_ij;
            SlipSystem(si).beta = beta_ij;
            SlipSystem(si).q = q_ij;
    end
%     save('FCC_SSSet.mat','SlipSystem');