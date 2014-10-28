function c = phys_const()
    r1 = 'CODATA 2010 (http://arxiv.org/abs/1203.5425)';
    
    c.c = 299792458;
    c.c_unit = 'm/s';
    c.c_err = 0;
    c.c_SI = c.c;
    c.c_unit_SI = c.c_unit;
    c.c_err_SI = 0;
    c.c_ref = r1;
    
    c.hbarc = 197.3269718;
    c.bharc_unit = 'Mev*fm';
    c.hbarc_err = 0.0000044;
    c.hbarc_ref = r1;
    
    c.mp = 938.272046;
    c.mp_unit = 'MeV';
    c.mp_err = 0.000021;
    c.mp_SI = 1.672621777e-27;
    c.mp_unit_SI = 'kg';
    c.mp_err_SI = 0.000000074e-27;
    c.mp_ref = r1;

    c.mn = 939.565379;
    c.mn_unit = 'MeV';
    c.mn_err = 0.000021;
    c.mn_SI = 1.674927351e-27;
    c.mn_unit_SI = 'kg';
    c.mn_err_SI = 0.000000074e-27;
    c.mn_ref = r1;
    
    c.md = 1875.612859;
    c.md_unit = 'MeV';
    c.md_err = 0.000041;
    c.md_SI = 3.34358348e-27;
    c.md_unit_SI = 'kg';
    c.md_err_SI = 0.00000015e-27;
    c.md_ref = r1;
    
    c.mt = 2808.921005;
    c.mt_unit = 'MeV';
    c.mt_err = 0.000062;
    c.mt_SI = 5.00735630e-27;
    c.mt_unit_SI = 'kg';
    c.mt_err_SI = 0.00000022e-27;
    c.mt_ref = r1;
    
    c.mh = 2808.391482;
    c.mh_unit = 'MeV';
    c.mh_err = 0.000062;
    c.mh_SI = 5.00641234e-27;
    c.mh_unit_SI = 'kg';
    c.mh_err_SI = 0.00000022e-27;
    c.mh_ref = r1;
    
    c.mHe = 3727.379240;
    c.mHe_unit = 'MeV';
    c.mHe_err = 0.000082;
    c.mHe_SI = 6.64465675e-27;
    c.mHe_unit_SI = 'kg';
    c.mHe_err_SI = 0.00000029e-27;
    c.mHe_ref = r1;
    
    c.mn_over_mp = 1.00137841917;
    c.mn_over_mp_unit = '';
    c.mn_over_mp_err = 0.00000000045;
    c.mn_over_mp_unit_SI = '';
    c.mn_over_mp_err_SI = c.mn_over_mp_err;
    c.mn_over_mp_ref = r1;
    
    c.md_over_mp = 1.99900750097;
    c.md_over_mp_unit = '';
    c.md_over_mp_err = 0.00000000018;
    c.md_over_mp_unit_SI = '';
    c.md_over_mp_err_SI = c.md_over_mp_err;
    c.md_over_mp_ref = r1;
    
    c.mt_over_mp = 2.9937170308;
    c.mt_over_mp_unit = '';
    c.mt_over_mp_err = 0.0000000025;
    c.mt_over_mp_unit_SI = '';
    c.mt_over_mp_err_SI = c.mt_over_mp_err;
    c.mt_over_mp_ref = r1;
    
    c.mh_over_mp = 2.9931526707;
    c.mh_over_mp_unit = '';
    c.mh_over_mp_err = 0.0000000025;
    c.mh_over_mp_unit_SI = '';
    c.mh_over_mp_err_SI = c.mh_over_mp_err;
    c.mh_over_mp_ref = r1;
    
    c.mHe_over_mp = 3.97259968933;
    c.mHe_over_mp_unit = '';
    c.mHe_over_mp_err = 0.00000000036;
    c.mHe_over_mp_unit_SI = '';
    c.mHe_over_mp_err_SI = c.mHe_over_mp_err;
    c.mHe_over_mp_ref = r1;
    
    c.u = 931.494061;
    c.u_unit = 'MeV';
    c.u_err = 0.000021;
    c.u_SI = 1.660538921e-27;
    c.u_unit_SI = 'kg';
    c.u_err_SI = 0.000000073e-27;
    c.u_ref = r1;
    
    c.me = 0.510998928;
    c.me_unit = 'MeV';
    c.me_err = 0.000000011;
    c.me_SI = 9.10938291e-31;
    c.me_unit_SI = 'kg';
    c.me_err_SI = 0.00000040e-31;
    c.me_ref = r1;
    
    c.atom_ionization_E_H_1 = 13.598434005136;
    c.atom_ionization_E_H_1_unit = 'eV';
    c.atom_ionization_E_H_1_err = 0.000000000012;
    c.atom_ionization_E_H_1_ref = 'NIST ASD (http://physics.nist.gov/asd)';
    
    c.a0 = 0.52917721092;
    c.a0_unit = 'Å';
    c.a0_err = 0.00000000017;
    c.a0_SI = c.a0*1e-10;
    c.a0_unit_SI = 'm';
    c.a0_err_SI = c.a0_err*1e-10;
    c.a0_ref = r1;
end
