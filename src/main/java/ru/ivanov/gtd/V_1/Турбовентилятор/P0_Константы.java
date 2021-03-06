package ru.ivanov.gtd.V_1.Турбовентилятор;

/**
 * @author Sergey Ivanov
 * created on 05.10.2021
 */
public class P0_Константы {

    public static final double otn_d_vt_1 = 0.35;
    public static final double sigma_vx = 1;
    public static final double sigma_na = 1;

    public static final double k_vozdyxa = 1.4;
    public static final double k_gaza = 1.33;
    public static final double R_vozdyxa = 287.3;
    public static final double R_gaza = 288.3;
    public static final double Nu_v = 0.85;
    public static final double Nu_kvd = 0.84;
    public static final double Nu_v_1 = 0.88;
    public static final double Nu_tvd = 0.88;
    public static final double Nu_tv$ = 0.9;

    public static final double s_v = 0.0404;
    public static final double s_g = 0.0396;
    public static final double Kg = 0.93;
    public static final double u_v_1 = 360;
    public static final double alfa_1_vblx = 80;
    public static final double sigma_pereh = 0.99;

    //компрессор выбирается из промежутка для получения скорости 180-220
    public static final double lambda_1_а = 0.53;

    //stepen reaktivnosti stypeni
    public static final double Rou_st_1 = 0.5;
    public static final double Rou_st_2 = 0.53;
    public static final double Rou_st_3 = 0.56;

    //koef koncevix poter
    public static final double Nu_konc_1 = 1;
    public static final double Nu_konc_2 = 0.99;
    public static final double Nu_konc_3 = 0.98;

    public static final double Nu_rk = 0.93;

    //выбирается в пределах
    public static final double lambda_2_tv = 0.48;
    //в градусах
    public static final double alfa_2_t = 87;

    //основные функции гтд
    public static final double lambda_vx = 0.65;
    public static final double q_Lambda_vx = 0.8543;

    public static final double lambda_vx_1 = 0.45;
    public static final double q_Lambda_vx_1 = 0.6515;
}
