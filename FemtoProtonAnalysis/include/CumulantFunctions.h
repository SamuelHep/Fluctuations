#ifndef CUMULANT_FUNCTIONS_H
#define CUMULANT_FUNCTIONS_H

long double q_ab( long double m_rs, long double m_r, long double m_s );
long double q_abc( long double m_rst, long double m_rs, long double m_rt, long double m_st, long double m_r, long double m_s, long double m_t );
long double q_abcd( long double m_rstu,
		    long double m_rst, long double m_rsu, long double m_rtu, long double m_stu,
		    long double m_rs, long double m_rt, long double m_ru, long double m_st, long double m_su, long double m_tu,
		    long double m_r, long double m_s, long double m_t, long double m_u
		    ,bool coutON=false);
long double q_abcde( long double m_rstuv,
		     long double m_rstu, long double m_rstv, long double m_rsuv, long double m_rtuv, long double m_stuv,
		     long double m_rst,long double m_rsu,long double m_rtu,long double m_stu,long double m_rsv,
		     long double m_rtv,long double m_stv,long double m_ruv,long double m_suv,long double m_tuv,
		     long double m_rs,long double m_rt,long double m_st,long double m_ru,long double m_su,
		     long double m_tu,long double m_rv,long double m_sv,long double m_tv,long double m_uv,
		     long double m_r, long double m_s, long double m_t, long double m_u, long double m_v
		     ,bool coutON=false);

#endif
