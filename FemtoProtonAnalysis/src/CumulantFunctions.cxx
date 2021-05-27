#include "CumulantFunctions.h"

#include <iostream>

using namespace std;

long double q_ab( long double m_rs, long double m_r, long double m_s )
{
  return m_rs - m_r*m_s;
} 

long double q_abc( long double m_rst, long double m_rs, long double m_rt, long double m_st, long double m_r, long double m_s, long double m_t )
{
  return m_rst - q_ab( m_rs, m_r, m_s )*m_t - q_ab( m_rt, m_r, m_t )*m_s - q_ab( m_rt, m_r, m_t )*m_s - m_r*m_s*m_t;
}

long double q_abcd( long double m_rstu,
		    long double m_rst, long double m_rsu, long double m_rtu, long double m_stu,
		    long double m_rs, long double m_rt, long double m_ru, long double m_st, long double m_su, long double m_tu,
		    long double m_r, long double m_s, long double m_t, long double m_u
		    ,bool coutON)
{
  long double q_rst = q_abc( m_rst, m_rs, m_rt, m_st, m_r, m_s, m_t );
  long double q_rsu = q_abc( m_rsu, m_rs, m_ru, m_su, m_r, m_s, m_u );
  long double q_rtu = q_abc( m_rtu, m_rt, m_ru, m_tu, m_r, m_t, m_u );
  long double q_stu = q_abc( m_stu, m_st, m_su, m_tu, m_s, m_t, m_u );

  long double q_rs = q_ab( m_rs, m_r, m_s);
  long double q_rt = q_ab( m_rt, m_r, m_t);
  long double q_ru = q_ab( m_ru, m_r, m_u);
  long double q_st = q_ab( m_st, m_s, m_t);
  long double q_su = q_ab( m_su, m_s, m_u);
  long double q_tu = q_ab( m_tu, m_t, m_u);

  long double q_r = m_r;
  long double q_s = m_s;
  long double q_t = m_t;
  long double q_u = m_u;

  if ( coutON )
    {

      cout << "m_rtu = " << m_rtu << endl;
      cout << "m_rt = " << m_rt << endl;
      cout << "m_ru = " << m_ru << endl;
      cout << "m_tu = " << m_tu << endl;
      cout << "m_r = " << m_r << endl;
      cout << "m_t = " << m_t << endl;
      cout << "m_u = " << m_u << endl;


	cout << "moment = " << m_rstu  << endl;
      cout << "term0 = " << q_rst*q_u << endl;
      cout << "term1 = " << q_rsu*q_t  << endl;
      cout << "term2 = " << q_rtu*q_s << " " << q_rtu << " " << q_s << endl;
      cout << "term3 = " << q_stu*q_r << endl;
      cout << "term4 = " << q_rs*q_tu  << endl;
      cout << "term5 = " << q_rt*q_su  << endl;
      cout << "term6 = " << q_ru*q_st  << endl;
      cout << "term7 = " << q_rs*q_t*q_u  << endl;
      cout << "term8 = " << q_rt*q_s*q_u << endl;
      cout << "term9 = " << q_ru*q_s*q_t << endl;
      cout << "term10 = " << q_st*q_r*q_u  << endl;
      cout << "term11 = " << q_su*q_r*q_t << endl;
      cout << "term12 = " << q_tu*q_r*q_s << endl;
      cout << "term13 = " << q_r*q_s*q_t*q_u << endl;

    }


  long double q_rstu = m_rstu 
    - q_rst*q_u - q_rsu*q_t - q_rtu*q_s - q_stu*q_r
    - q_rs*q_tu - q_rt*q_su  - q_ru*q_st 
    - q_rs*q_t*q_u - q_rt*q_s*q_u - q_ru*q_s*q_t - q_st*q_r*q_u - q_su*q_r*q_t - q_tu*q_r*q_s
    -q_r*q_s*q_t*q_u;

  return q_rstu;

}


// Cumulants to Moments
// m_rstu = q_rstu + q_rst*qu[4] + q_rs*q_tu[3] + q_rs*q_t*_qu[6] + q_r*q_s*q_t*q_u 
// 4th Cumulant
// q_rstu = m_rstu - q_rst*qu[4] - q_rs*q_tu[3] - q_rs*q_t*q_u[6] - q_r*q_s*q_t*q_u 

// 5th Moment
// m_rstuv = q_rstuv + q_rstu*q_v[5] + q_rst*q_uv[10] + q_rst*q_u*q_v[10] + q_rs*q_t*q_u*q_v[10] +  q_rs*q_tu*q_v[15] + q_r*q_s*q_t*q_u*q_v
// 5th Cumulant
// q_rstuv = m_rstuv - q_rstu*q_v[5] - q_rst*q_uv[10] - q_rst*q_u*q_v[10] - q_rs*q_t*q_u*q_v[10] -  q_rs*q_tu*q_v[15] - q_r*q_s*q_t*q_u*q_v

long double q_abcde( long double m_rstuv,
		     long double m_rstu, long double m_rstv, long double m_rsuv, long double m_rtuv, long double m_stuv,
		     long double m_rst,long double m_rsu,long double m_rtu,long double m_stu,long double m_rsv,
		     long double m_rtv,long double m_stv,long double m_ruv,long double m_suv,long double m_tuv,
		     long double m_rs,long double m_rt,long double m_st,long double m_ru,long double m_su,
		     long double m_tu,long double m_rv,long double m_sv,long double m_tv,long double m_uv,
		     long double m_r, long double m_s, long double m_t, long double m_u, long double m_v
		     ,bool coutON )

{

  long double q_rstu = q_abcd( m_rstu, m_rst, m_rsu, m_rtu, m_stu,m_rs, m_rt, m_ru, m_st, m_su, m_tu, m_r, m_s, m_t, m_u );
  long double q_rstv = q_abcd( m_rstv, m_rst, m_rsv, m_rtv, m_stv,m_rs, m_rt, m_rv, m_st, m_sv, m_tv, m_r, m_s, m_t, m_v );
  long double q_rsuv = q_abcd( m_rsuv, m_rsu, m_rsv, m_ruv, m_suv,m_rs, m_ru, m_rv, m_su, m_sv, m_uv, m_r, m_s, m_u, m_v );
  long double q_rtuv = q_abcd( m_rtuv, m_rtu, m_rtv, m_ruv, m_tuv,m_rt, m_ru, m_rv, m_tu, m_tv, m_uv, m_r, m_t, m_u, m_v );
  long double q_stuv = q_abcd( m_stuv, m_stu, m_stv, m_suv, m_tuv,m_st, m_su, m_sv, m_tu, m_tv, m_uv, m_s, m_t, m_u, m_v );

  long double q_rst = q_abc( m_rst, m_rs, m_rt, m_st, m_r, m_s, m_t );
  long double q_rsu = q_abc( m_rsu, m_rs, m_ru, m_su, m_r, m_s, m_u );
  long double q_rtu = q_abc( m_rtu, m_rt, m_ru, m_tu, m_r, m_t, m_u );
  long double q_stu = q_abc( m_stu, m_st, m_su, m_tu, m_s, m_t, m_u );
  long double q_rsv = q_abc( m_rsv, m_rs, m_rv, m_sv, m_r, m_s, m_v );
  
  long double q_rtv = q_abc( m_rtv, m_rt, m_rv, m_tv, m_r, m_t, m_v );  
  long double q_stv = q_abc( m_stv, m_st, m_sv, m_tv, m_s, m_t, m_v );
  long double q_ruv = q_abc( m_ruv, m_ru, m_rv, m_uv, m_r, m_u, m_v );
  long double q_suv = q_abc( m_suv, m_su, m_sv, m_uv, m_s, m_u, m_v );
  long double q_tuv = q_abc( m_tuv, m_tu, m_tv, m_uv, m_t, m_u, m_v );

  long double q_rs = q_ab( m_rs, m_r, m_s);
  long double q_rt = q_ab( m_rt, m_r, m_t);
  long double q_st = q_ab( m_st, m_s, m_t);
  long double q_ru = q_ab( m_ru, m_r, m_u);
  long double q_su = q_ab( m_su, m_s, m_u);

  long double q_tu = q_ab( m_tu, m_t, m_u);
  long double q_rv = q_ab( m_rv, m_r, m_v);
  long double q_sv = q_ab( m_sv, m_s, m_v);
  long double q_tv = q_ab( m_tv, m_t, m_v);
  long double q_uv = q_ab( m_uv, m_u, m_v);
  
  long double q_r = m_r;
  long double q_s = m_s;
  long double q_t = m_t;
  long double q_u = m_u;
  long double q_v = m_v;

  long double term1 = m_rstuv;
  long double term2 = q_rstu*q_v + q_rstv*q_u + q_rsuv*q_t + q_rtuv*q_s + q_stuv*q_r;
  long double term3 = q_st*q_ruv + q_su*q_rtv + q_tv*q_rsu + q_sv*q_rtu + q_rv*q_stu + q_tu*q_rsv + q_rs*q_tuv + q_rt*q_suv + q_ru*q_stv + q_uv*q_rst;
  long double term4 = q_s*q_t*q_ruv + q_s*q_u*q_rtv + q_t*q_v*q_rsu + q_s*q_v*q_rtu + q_r*q_v*q_stu + q_t*q_u*q_rsv + q_r*q_s*q_tuv + q_r*q_t*q_suv + q_r*q_u*q_stv + q_u*q_v*q_rst;
  long double term5 = q_st*q_r*q_u*q_v + q_su*q_r*q_t*q_v + q_tv*q_r*q_s*q_u + q_sv*q_r*q_t*q_u + q_rv*q_s*q_t*q_u + q_tu*q_r*q_s*q_v + q_rs*q_t*q_u*q_v + q_rt*q_s*q_u*q_v + q_ru*q_s*q_t*q_v + q_uv*q_r*q_s*q_t;
  long double term6 = q_uv*q_st*q_r + q_rt*q_su*q_v + q_rt*q_sv*q_u + q_sv*q_tu*q_r + q_uv*q_rs*q_t + q_uv*q_rt*q_s + q_rv*q_st*q_u + q_tv*q_rs*q_u + q_tv*q_su*q_r + q_ru*q_sv*q_t + q_ru*q_st*q_v + q_rv*q_tu*q_s + q_ru*q_tv*q_s + q_rv*q_su*q_t + q_rs*q_tu*q_v;
  long double term7 = q_r*q_s*q_t*q_u*q_v;


  if ( coutON )
    {
      cout << "TERM 2 ... " << endl;
      cout << "q_rstu " << q_rstu << endl;
      cout << "q_rstv " << q_rstv << endl;
      cout << "q_rsuv " << q_rsuv << endl;
      cout << "q_rtuv " << q_rtuv << endl;
      cout << "q_stuv " << q_stuv << endl;

      cout << "moment = " << m_rstuv << endl;
      cout << "TERM 1 = " << term1 << endl;
      cout << "TERM 2 = " << term2 << endl;
      cout << "TERM 3 = " << term3 << endl;
      cout << "TERM 4 = " << term4 << endl;
      cout << "TERM 5 = " << term5 << endl;
      cout << "TERM 6 = " << term6 << endl;
      cout << "TERM 7 = " << term7 << endl;
    }

  return term1 - term2 - term3 - term4 - term5 - term6 - term7;

  

}



