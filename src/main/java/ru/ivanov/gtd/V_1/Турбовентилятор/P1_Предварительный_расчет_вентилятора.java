package ru.ivanov.gtd.V_1.Турбовентилятор;

import static ru.ivanov.gtd.V_1.Util.*;
import static ru.ivanov.gtd.V_1.Турбовентилятор.P00_Исходные_данные.*;
import static ru.ivanov.gtd.V_1.Турбовентилятор.P0_Константы.*;

/**
 * @author Sergey Ivanov
 * created on 05.10.2021
 */
public class P1_Предварительный_расчет_вентилятора {
    public double P_vx$;
    public double D_v_1;
    public double D_sr_1;
    public double D_vt_1;

    public double otn_F;
    public double L_v;
    public double T_k_vx$, T_v_vblx$;
    public double n; // показатель политропного процесса

    public double F_v_vblx;
    public double otn_d_vblx;

    public double D_vblx_nar;
    public double D_vblx_sr;
    public double D_vblx_vt;
    public double sum_otn_L_v;

    public double otn_L_v_1;
    public double otn_L_v_2;
    public double otn_L_v_3;

    public double L_v_1;
    public double L_v_2;
    public double L_v_3;

    public double sum_L_v;
    public double u_sr_1;


    private void calculate_P_vx$() {
        this.P_vx$ = P_n$ * sigma_vx * sigma_na;
    }

    private void calculate_D_v_1() {
        double numerator = G_v_sum * Math.sqrt(T_vx$);
        double denominator = Math.PI * s_v * P_vx$ * (1 - otn_d_vt_1 * otn_d_vt_1) * q_Lambda_vx * Kg;
        this.D_v_1 = 2 * Math.sqrt(numerator / denominator);
        System.out.println("D_v_1 = " + D_v_1);
    }

    private void calculate_D_sr_1() {
        double arg = (1 + otn_d_vt_1 * otn_d_vt_1) / 2;
        this.D_sr_1 = D_v_1 * Math.sqrt(arg);
        System.out.println("D_sr_1 = " + D_sr_1);
    }

    private void calculate_D_vt_1() {
        this.D_vt_1 = otn_d_vt * D_v_1;
        System.out.println("D_vt_1 = " + D_vt_1);
    }

    private void calculate_otn_F() {
        this.L_v = (k_vozdyxa / (k_vozdyxa - 1)) * R_vozdyxa * T_n$ * (Math.pow(P00_Исходные_данные.Pi_v$, (k_vozdyxa - 1) / k_vozdyxa) - 1) * (1 / Nu_v);
        this.T_k_vx$ = T_vx$ + (L_v / ((k_vozdyxa / (k_vozdyxa - 1)) * R_vozdyxa));
        this.T_v_vblx$ = T_k_vx$;
        double arg = Math.log(P00_Исходные_данные.Pi_v$) / Math.log(T_k_vx$ / T_vx$);
        this.n = arg / (arg - 1);
        // q_Lambda_k_vx / q_Lambda_vx = 1
        this.otn_F = Math.pow(P00_Исходные_данные.Pi_v$, (n + 1) / (2 * n));
        System.out.println("L_v = " + L_v);
        System.out.println("T_k_vx$ = " + T_k_vx$);
        System.out.println("T_v_vblx$ = " + T_v_vblx$);
        System.out.println("n = " + n);
        System.out.println("otn_F = " + otn_F);
    }

    private void calculate_F_v_vblx() {
        double numerator = G_v_sum * Math.sqrt(T_vx$);
        double denominator = s_v * P_vx$ * q_Lambda_vx * Kg * otn_F;
        this.F_v_vblx = numerator / denominator;
        System.out.println("F_v_vblx = " + F_v_vblx);
    }

    private void calculate_otn_d_vblx() {
        this.otn_d_vblx = Math.sqrt((otn_F + otn_d_vt_1 * otn_d_vt_1 - 1) / (otn_F));
        System.out.println("otn_d_vblx = " + otn_d_vblx);
    }

    private void calculate_D_vblx_sr() {
        this.D_vblx_nar = D_v_1;
        this.D_vblx_sr = D_v_1 * Math.sqrt((1 + otn_d_vblx * otn_d_vblx) / 2);
        System.out.println("D_vblx_nar = " + D_vblx_nar);
        System.out.println("D_vblx_sr = " + D_vblx_sr);
    }

    private void calculate_D_vblx_vt() {
        this.D_vblx_vt = D_vblx_nar * otn_d_vblx;
        System.out.println("D_vblx_vt = " + D_vblx_vt);
    }

    private void calculate_sum_otn_L_v() {
        this.sum_otn_L_v = L_v / (u_v_1 * u_v_1);
        System.out.println("sum_otn_L_v = " + sum_otn_L_v);
    }

    private void PODBOR_STYPENEY() {
        //алгоритм подбора нужен сюда
        this.otn_L_v_1 = 0.29;
        this.otn_L_v_2 = 0.346;
        this.otn_L_v_3 = 0.4;
        System.out.println("otn_L_v_1 = " + otn_L_v_1);
        System.out.println("otn_L_v_2 = " + otn_L_v_2);
        System.out.println("otn_L_v_3 = " + otn_L_v_3);
    }

    private void calculate_L_v() {
        this.L_v_1 = otn_L_v_1 * u_v_1 * u_v_1;
        this.L_v_2 = otn_L_v_2 * u_v_1 * u_v_1;
        this.L_v_3 = otn_L_v_3 * u_v_1 * u_v_1;
        System.out.println("L_v_1 = " + L_v_1);
        System.out.println("L_v_2 = " + L_v_2);
        System.out.println("L_v_3 = " + L_v_3);
    }

    private void calculate_sum_L_v() {
        this.sum_L_v = L_v_1 + L_v_2 + L_v_3;
        System.out.println("sum_L_v = " + sum_L_v);
    }

    private void calculate_u_sr_1() {
        this.u_sr_1 = u_v_1 * Math.sqrt((1 + otn_d_vt_1 * otn_d_vt_1) / 2);
        System.out.println("u_sr_1 = " + u_sr_1);
    }


    public double otn_L£_v_1;
    public double otn_L£_v_2;
    public double otn_L£_v_3;
    public double D_sr_2;
    public double D_sr_3;
    public double u_sr_2;
    public double u_sr_3;

    //£
    private void koef_nagryzki_stypenei() {
        this.otn_L£_v_1 = L_v_1 / (u_sr_1 * u_sr_1);
        this.D_sr_2 = (D_sr_1 + D_vblx_sr) / 2;
        this.u_sr_2 = u_sr_1 * (D_sr_2 / D_sr_1);
        this.otn_L£_v_2 = L_v_2 / (u_sr_2 * u_sr_2);
        this.D_sr_3 = (D_sr_2 + D_vblx_sr) / 2;
        this.u_sr_3 = u_sr_2 * (D_sr_3 / D_sr_2);
        this.otn_L£_v_3 = L_v_3 / (u_sr_3 * u_sr_3);
        System.out.println("otn_L£_v_1 = " + otn_L£_v_1);
        System.out.println("otn_L£_v_2 = " + otn_L£_v_2);
        System.out.println("otn_L£_v_3 = " + otn_L£_v_3);
        System.out.println("D_sr_2 = " + D_sr_2);
        System.out.println("D_sr_3 = " + D_sr_3);
        System.out.println("u_sr_2 = " + u_sr_2);
        System.out.println("u_sr_3 = " + u_sr_3);
    }

    /////////////////////////////////ТАБЛИЦА 1 (ПОСТРОЕНИЕ)

    public double N_v;

    public void calculate_N_v() {
        this.N_v = (u_sr_1 * 60) / (Math.PI * D_sr_1);
        System.out.println("N_v = " + N_v);
    }

    public double F_vblx_1;
    public double G_v_1;
    public double P_k_vx$;


    public void calculate_F_vblx_1() {
        this.G_v_1 = (1 / (1 + m)) * G_v_sum;
        this.P_k_vx$ = P00_Исходные_данные.Pi_v$ * P_n$ * ((sigma_vx * sigma_na) / sigma_pereh);
        double numerator = G_v_1 * Math.sqrt(T_v_vblx$) * sigma_na * sigma_pereh;
        double denominator = s_v * P_k_vx$ * q_Lambda_vx_1 * Kg * Math.sin(fromGradToRad(alfa_1_vblx));
        System.out.println(fromGradToRad(alfa_1_vblx));
        this.F_vblx_1 = numerator / denominator;
        System.out.println("G_v_1 = " + G_v_1);
        System.out.println("P_k_vx$ = " + P_k_vx$);
        System.out.println("F_vblx_1 = " + F_vblx_1);
    }

    public double F_vblx_2;

    public void calculate_F_vblx_2() {
        this.F_vblx_2 = F_v_vblx - F_vblx_1;
        System.out.println("F_vblx_2 = " + F_vblx_2);
    }

    public double D_razd;

    public void calculate_D_razd() {
        this.D_razd = Math.sqrt((D_vblx_nar * D_vblx_nar) - ((4 / Math.PI) * F_vblx_2));
        System.out.println("D_razd = " + D_razd);
    }

    public double T_3_1$;
    public double T_3_2$;
    public double T_3_3$;
    public double T_2_1$;
    public double T_2_2$;
    public double T_2_3$;
    public double T_1_1$;
    public double T_1_2$;
    public double T_1_3$;
    public double T_vblx$;

    private void calculate_izontrop_potok() {
        System.out.println("--------------");
        this.T_1_1$ = 288.0;
        this.T_1_2$ = T_1_1$ + (L_v_1 / (k_vozdyxa / (k_vozdyxa - 1) * R_vozdyxa));
        this.T_2_1$ = T_1_2$;
        this.T_3_1$ = T_1_2$;

        this.T_1_3$ = T_1_2$ + (L_v_2 / (k_vozdyxa / (k_vozdyxa - 1) * R_vozdyxa));
        this.T_2_2$ = T_1_3$;
        this.T_3_2$ = T_1_3$;

        this.T_2_3$ = T_1_3$ + (L_v_3 / (k_vozdyxa / (k_vozdyxa - 1) * R_vozdyxa));
        this.T_3_3$ = T_2_3$;
        this.T_vblx$ = T_3_3$;
        System.out.println("T_1_1$ = " + T_1_1$);
        System.out.println("T_1_2$ = " + T_1_2$);
        System.out.println("T_1_3$ = " + T_1_3$);

        System.out.println("T_2_1$ = " + T_2_1$);
        System.out.println("T_2_2$ = " + T_2_2$);
        System.out.println("T_2_3$ = " + T_2_3$);

        System.out.println("T_3_1$ = " + T_3_1$);
        System.out.println("T_3_2$ = " + T_3_2$);
        System.out.println("T_3_3$ = " + T_3_3$);
        System.out.println("T_vblx$ = " + T_vblx$);
        System.out.println("--------------");
    }

    public double Pi_v_1$;
    public double Pi_v_2$;
    public double Pi_v_3$;

    public void calculate_pi_v() {
        this.Pi_v_1$ = Math.pow((L_v_1 * Nu_v_1) / ((k_vozdyxa / (k_vozdyxa - 1)) * R_vozdyxa * T_1_1$) + 1, (k_vozdyxa / (k_vozdyxa - 1)));
        this.Pi_v_2$ = Math.pow((L_v_2 * Nu_v_1) / ((k_vozdyxa / (k_vozdyxa - 1)) * R_vozdyxa * T_1_2$) + 1, (k_vozdyxa / (k_vozdyxa - 1)));
        this.Pi_v_3$ = Math.pow((L_v_3 * Nu_v_1) / ((k_vozdyxa / (k_vozdyxa - 1)) * R_vozdyxa * T_1_3$) + 1, (k_vozdyxa / (k_vozdyxa - 1)));
        System.out.println("Pi_v_1$ = " + Pi_v_1$);
        System.out.println("Pi_v_2$ = " + Pi_v_2$);
        System.out.println("Pi_v_2$ = " + Pi_v_3$);
        double Pi_v$ = Pi_v_1$ * Pi_v_2$ * Pi_v_3$;
        System.out.println("Pi_v$ как произведение ПиВ на каждой ступени (для проверки) = " + Pi_v$);
        System.out.println("Pi_v$ изначальное = " + P00_Исходные_данные.Pi_v$);
        System.out.printf("Разница ПиВ = %.2f", (Pi_v$ - P00_Исходные_данные.Pi_v$) / P00_Исходные_данные.Pi_v$ * 100);
    }

    public double P_1_1$;
    public double P_1_2$;
    public double P_1_3$;
    public double P_3_1$;
    public double P_3_2$;
    public double P_3_3$;
    public double P_v_vblx$;

    public void polnoe_davlenie_na_vxode_v_stypeni() {
        this.P_1_1$ = P_n$ * sigma_na;
        this.P_3_1$ = P_1_1$ * Pi_v_1$;
        this.P_1_2$ = P_3_1$;
        this.P_3_2$ = P_1_2$ * Pi_v_2$;
        this.P_1_3$ = P_3_2$;
        this.P_3_3$ = P_1_3$ * Pi_v_3$;
        this.P_v_vblx$ = P_3_3$;
        System.out.println("\n-----------------------");
        System.out.println("P_1_1$ = " + P_1_1$);
        System.out.println("P_1_2$ = " + P_1_2$);
        System.out.println("P_1_3$ = " + P_1_3$);
        System.out.println("P_3_1$ = " + P_3_1$);
        System.out.println("P_3_2$ = " + P_3_2$);
        System.out.println("P_3_3$ = " + P_3_3$);
        System.out.println("P_v_vblx$ = " + P_v_vblx$);
        System.out.println("-----------------------");
    }

    public double alfa_1_kr_1;
    public double alfa_1_kr_2;
    public double alfa_1_kr_3;
    public double alfa_1_kr_vblx;

    public void kriticheskie_skorosti_na_vxode() {
        this.alfa_1_kr_1 = Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_1_1$);
        this.alfa_1_kr_2 = Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_1_2$);
        this.alfa_1_kr_3 = Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_1_3$);
        this.alfa_1_kr_vblx = Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_vblx$);
        System.out.println("\n-----------------------");
        System.out.println("alfa_1_kr_1 = " + alfa_1_kr_1);
        System.out.println("alfa_1_kr_2 = " + alfa_1_kr_2);
        System.out.println("alfa_1_kr_3 = " + alfa_1_kr_3);
        System.out.println("alfa_1_kr_vblx = " + alfa_1_kr_vblx);
        System.out.println("-----------------------");
    }

    public double c_1_a_1;
    public double c_1;

    public void skorost_na_vxode_v_ventilyator() {
        this.c_1_a_1 = c_1 = lambda_vx * alfa_1_kr_1;
        System.out.println("c_1_a_1 = " + c_1_a_1);
    }

    //kompressor dlya skorostey
    public double alfa_1_kr_kompressor;
    public double c_1_alfa_kompressor;

    public void osevaya_skorost_kompressora() {
        System.out.println("-----------KOMPRESSOR-----------");
        this.alfa_1_kr_kompressor = Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_k_vx$);
        this.c_1_alfa_kompressor = lambda_1_а * alfa_1_kr_kompressor;
        System.out.println("alfa_1_kr_kompressor = " + alfa_1_kr_kompressor);
        System.out.println("c_1_alfa_kompressor = " + c_1_alfa_kompressor);
        System.out.println("-----------KOMPRESSOR-----------");
    }

    public double c_1_a_2;
    public double c_1_a_3;
    public double c_1_a_vblx;

    public void skorosti_na_stypenyax() {
        //распределение от c_1_a_1 до c_1_a_vblx сделать
        this.c_1_a_2 = 201.5;
        this.c_1_a_3 = 199.8;
        this.c_1_a_vblx = c_1_alfa_kompressor;
        System.out.println("\n-----------------------");
        System.out.println("c_1_a_2 = " + c_1_a_2);
        System.out.println("c_1_a_3 = " + c_1_a_3);
        System.out.println("c_1_a_vblx = " + c_1_a_vblx);
        System.out.println("-----------------------");
    }

    public double otn_c_1_alfa_1;
    public double otn_c_1_alfa_2;
    public double otn_c_1_alfa_3;

    public void koef_rasxoda_na_srednem_diametre_kolesa() {
        this.otn_c_1_alfa_1 = c_1_a_1 / u_sr_1;
        this.otn_c_1_alfa_2 = c_1_a_2 / u_sr_2;
        this.otn_c_1_alfa_3 = c_1_a_3 / u_sr_3;
        System.out.println("\n-----------------------");
        System.out.println("otn_c_1_alfa_1 = " + otn_c_1_alfa_1);
        System.out.println("otn_c_1_alfa_2 = " + otn_c_1_alfa_2);
        System.out.println("otn_c_1_alfa_3 = " + otn_c_1_alfa_3);
        System.out.println("-----------------------");
    }

    public double otn_L_k_u_1;
    public double otn_L_k_u_2;
    public double otn_L_k_u_3;

    public void koef_napora_na_kolese() {
        this.otn_L_k_u_1 = (L_v_1 * Nu_konc_1 * Nu_konc_1) / (u_sr_1 * u_sr_1);
        this.otn_L_k_u_2 = (L_v_2 * Nu_konc_2 * Nu_konc_2) / (u_sr_2 * u_sr_2);
        this.otn_L_k_u_3 = (L_v_3 * Nu_konc_3 * Nu_konc_3) / (u_sr_3 * u_sr_3);
        System.out.println("\n-----------------------");
        System.out.println("otn_L_k_u_1 = " + otn_L_k_u_1);
        System.out.println("otn_L_k_u_2 = " + otn_L_k_u_2);
        System.out.println("otn_L_k_u_3 = " + otn_L_k_u_3);
        System.out.println("-----------------------");
    }

    //ygol vxoda vozdyxa v koleso
    public double katangens_alfa_1_1;
    public double alfa_1_1;

    public void ygol_cxods_vozdyxa_v_koleso_na_srednem_radiyse() {
        this.katangens_alfa_1_1 = (2 * (1 - Rou_st_1) - otn_L_k_u_1) / (2 * otn_c_1_alfa_1);
        this.alfa_1_1 = arcctg(katangens_alfa_1_1);
        System.out.println("\n-----------------------");
        System.out.println("katangens_alfa_1_1 = " + katangens_alfa_1_1);
        System.out.println("alfa_1_1 = " + alfa_1_1);
        System.out.println("-----------------------");
    }

    //prived skorost
    public double lamda_1_1;

    public void prived_skorost_na_vxode() {
        this.lamda_1_1 = c_1_a_1 / (alfa_1_kr_1 * Math.sin(fromGradToRad(alfa_1_1)));
        System.out.println("lamda_1_1 = " + lamda_1_1);
    }

    public double M_w_1;
    public double lambda_1_u;
    public double tau_lamda_1_1 = 0.9229;

    // тут альфы критические и альфы СОВСЕМ ДРУГИЕ
    public void chislo_maxa_na_vxode_po_otn_skorosti() {

        //ne alfa - drygoy koef
        double alfa_1_kr_1 = 1;
        //ne alfa - drygoy koef
        double alfa_1_1 = Math.sqrt(((k_vozdyxa + 1) / 2) * tau_lamda_1_1);

        this.lambda_1_u = u_sr_1 / this.alfa_1_kr_1;

        this.M_w_1 = (alfa_1_kr_1 / alfa_1_1) * Math.sqrt(lamda_1_1 * lamda_1_1 + lambda_1_u * lambda_1_u - 2 * lamda_1_1 * lambda_1_u * Math.cos(fromGradToRad(this.alfa_1_1)));
        System.out.println("tau_lamda_1_1 = " + tau_lamda_1_1);
        System.out.println("alfa_1_1 = " + alfa_1_1);
        System.out.println("(alfa_1_kr_1 / alfa_1_1) = " + (alfa_1_kr_1 / alfa_1_1));
        System.out.println("lambda_1_u = " + lambda_1_u);
        System.out.println("M_w_1 = " + M_w_1);
        System.out.println("-----------------------");
    }

    public double otn_G_k_1;
    public double q_lamda_1_1 = 0.8778;

    public void koef_proizvoditelnosti_pervoi_stypeni() {
        this.otn_G_k_1 = (1 - otn_d_vt_1 * otn_d_vt_1) * q_lamda_1_1 * Math.sin(fromGradToRad(alfa_1_1));
        System.out.println("otn_G_k_1 = " + otn_G_k_1);
    }

    public double c_1_u_1;
    public double c_1_u_2;
    public double c_1_u_3;
    public double Rou_k_1;
    public double Rou_k_2;
    public double Rou_k_3;

    public void okryzhnaya_sostavl_absolut_skorosti() {
        this.Rou_k_1 = Rou_st_1;
        this.Rou_k_2 = Rou_st_2;
        this.Rou_k_3 = Rou_st_3;
        this.c_1_u_1 = u_sr_1 * ((1 - Rou_k_1) - otn_L_k_u_1 / 2);
        this.c_1_u_2 = u_sr_2 * ((1 - Rou_k_2) - otn_L_k_u_2 / 2);
        this.c_1_u_3 = u_sr_3 * ((1 - Rou_k_3) - otn_L_k_u_3 / 2);
        System.out.println("\n-----------------------");
        System.out.println("Rou_k_1 = " + Rou_k_1);
        System.out.println("Rou_k_2 = " + Rou_k_2);
        System.out.println("Rou_k_3 = " + Rou_k_3);
        System.out.println("c_1_u_1 = " + c_1_u_1);
        System.out.println("c_1_u_2 = " + c_1_u_2);
        System.out.println("c_1_u_3 = " + c_1_u_3);
        System.out.println("-----------------------");
    }

    public double c_1_1;
    public double c_1_2;
    public double c_1_3;
    //tyt tozhe drygie lamda (lamda_1_1 dolzhno bit);
    double lamda_1_1_f;
    double lamda_1_2_f;
    double lamda_1_3_f;
    double lamda_1_vblx_f;

    //DOBAVIL f postsyfix
    public void absolut_i_prived_skorosti_na_vxode_v_rab_kolesa() {
        this.c_1_1 = Math.sqrt(c_1_u_1 * c_1_u_1 + c_1_a_1 * c_1_a_1);
        this.c_1_2 = Math.sqrt(c_1_u_2 * c_1_u_2 + c_1_a_2 * c_1_a_2);
        this.c_1_3 = Math.sqrt(c_1_u_3 * c_1_u_3 + c_1_a_3 * c_1_a_3);
        this.lamda_1_1_f = c_1_1 / Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_1_1$);
        this.lamda_1_2_f = c_1_2 / Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_1_2$);
        this.lamda_1_3_f = c_1_3 / Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_1_3$);
        this.lamda_1_vblx_f = c_1_3 / Math.sqrt(2 * k_vozdyxa / (k_vozdyxa + 1) * R_vozdyxa * T_vblx$);
        System.out.println("\n-----------------------");
        System.out.println("c_1_1 = " + c_1_1);
        System.out.println("c_1_2 = " + c_1_2);
        System.out.println("c_1_3 = " + c_1_3);
        System.out.println("lamda_1_1_f = " + lamda_1_1_f);
        System.out.println("lamda_1_2_f = " + lamda_1_2_f);
        System.out.println("lamda_1_3_f = " + lamda_1_3_f);
        System.out.println("lamda_1_vblx_f = " + lamda_1_vblx_f);
        System.out.println("-----------------------");
    }

    //tyt tozhe drygie lamda (alfa_1_1 dolzhno bit);
    double alfa_1_1_f;
    double alfa_1_2_f;
    double alfa_1_3_f;
    double alfa_v_vblx;

    public void ygol_vxoda_v_stypen_po_absolut_skorosti() {
        this.alfa_1_1_f = fromRadToGrad(Math.asin(c_1_a_1 / c_1_1));
        this.alfa_1_2_f = fromRadToGrad(Math.asin(c_1_a_2 / c_1_2));
        this.alfa_1_3_f = fromRadToGrad(Math.asin(c_1_a_3 / c_1_3));
        this.alfa_v_vblx = 90;
        System.out.println("\n-----------------------");
        System.out.println("alfa_1_1_f = " + alfa_1_1_f);
        System.out.println("alfa_1_2_f = " + alfa_1_2_f);
        System.out.println("alfa_1_3_f = " + alfa_1_3_f);
        System.out.println("alfa_v_vblx = " + alfa_v_vblx);
        System.out.println("-----------------------");
    }

    //tyt tozhe drygie lamda (q_lambda_1_1 dolzhno bit);
    double q_lamda_1_1_f = 0.8778;
    double q_lamda_1_2_f = 0.8459;
    double q_lamda_1_3_f = 0.8015;
    double q_lamda_1_vblx_f = 0.7623;
    //площаи на входе в ступени
    double F_1_1_f;
    double F_1_2_f;
    double F_1_3_f;
    double F_v_vblx_f;

    public void calculate_ploshad_proxodnogo_sechenia() {
        this.F_1_1_f = (G_v_sum * Math.sqrt(T_1_1$)) /
                (P_1_1$ * q_lamda_1_1_f * Math.sin(fromGradToRad(alfa_1_1_f)) * s_v * Kg);
        this.F_1_2_f = (G_v_sum * Math.sqrt(T_1_2$)) /
                (P_1_2$ * q_lamda_1_2_f * Math.sin(fromGradToRad(alfa_1_2_f)) * s_v * Kg);
        this.F_1_3_f = (G_v_sum * Math.sqrt(T_1_3$)) /
                (P_1_3$ * q_lamda_1_3_f * Math.sin(fromGradToRad(alfa_1_3_f)) * s_v * Kg);
        this.F_v_vblx_f = (G_v_sum * Math.sqrt(T_v_vblx$)) /
                (P_v_vblx$ * q_lamda_1_vblx_f * Math.sin(fromGradToRad(alfa_v_vblx)) * s_v * Kg);
        System.out.println("\n-----------------------");
        System.out.println("q_lamda_1_1_f = " + q_lamda_1_1_f);
        System.out.println("q_lamda_1_2_f = " + q_lamda_1_2_f);
        System.out.println("q_lamda_1_3_f = " + q_lamda_1_3_f);
        System.out.println("q_lamda_1_vblx_f = " + q_lamda_1_vblx_f);
        System.out.println("F_1_1_f = " + F_1_1_f);
        System.out.println("F_1_2_f = " + F_1_2_f);
        System.out.println("F_1_3_f = " + F_1_3_f);
        System.out.println("F_v_vblx_f = " + F_v_vblx_f);
        System.out.println("-----------------------");
    }

    public double D_v_1_1;
    public double D_v_1_2;
    public double D_v_1_3;
    public double D_v_vblx;

    public void naryzniy_diametr_rab_koles() {
        this.D_v_1_1 = Math.sqrt(D_sr_1 * D_sr_1 + (2 * F_1_1_f / Math.PI));
        this.D_v_1_2 = Math.sqrt(D_sr_2 * D_sr_2 + (2 * F_1_2_f / Math.PI));
        this.D_v_1_3 = Math.sqrt(D_sr_3 * D_sr_3 + (2 * F_1_3_f / Math.PI));
        this.D_v_vblx = Math.sqrt(D_vblx_sr * D_vblx_sr + (2 * F_v_vblx_f / Math.PI));

        System.out.println("\n-----------------------");
        System.out.println("D_v_1_1 = " + D_v_1_1);
        System.out.println("D_v_1_2 = " + D_v_1_2);
        System.out.println("D_v_1_3 = " + D_v_1_3);
        System.out.println("D_v_vblx = " + D_v_vblx);
        System.out.println("-----------------------");
    }

    public double D_vt_1_1;
    public double D_vt_1_2;
    public double D_vt_1_3;
    public double D_vt_vblx;

    public void vnytrenniy_diametr_rab_koles() {
        this.D_vt_1_1 = Math.sqrt(D_sr_1 * D_sr_1 - (2 * F_1_1_f / Math.PI));
        this.D_vt_1_2 = Math.sqrt(D_sr_2 * D_sr_2 - (2 * F_1_2_f / Math.PI));
        this.D_vt_1_3 = Math.sqrt(D_sr_3 * D_sr_3 - (2 * F_1_3_f / Math.PI));
        this.D_vt_vblx = Math.sqrt(D_vblx_sr * D_vblx_sr - (2 * F_v_vblx_f / Math.PI));

        System.out.println("\n-----------------------");
        System.out.println("D_vt_1_1 = " + D_vt_1_1);
        System.out.println("D_vt_1_2 = " + D_vt_1_2);
        System.out.println("D_vt_1_3 = " + D_vt_1_3);
        System.out.println("D_vt_vblx = " + D_vt_vblx);
        System.out.println("-----------------------");
    }

    //Dk = const значит D_v = D_v_vblx
    public double D_v;
    public double h_1_1;
    public double h_1_2;
    public double h_1_3;

    public void vblsota_rabochey_lopatki() {
        this.D_v = D_v_vblx;
        this.h_1_1 = (D_v - D_vt_1_1) / 2;
        this.h_1_2 = (D_v - D_vt_1_2) / 2;
        this.h_1_3 = (D_v - D_vt_1_3) / 2;

        System.out.println("\n-----------------------");
        System.out.println("D_v = " + D_v);
        System.out.println("h_1_1 = " + h_1_1);
        System.out.println("h_1_2 = " + h_1_2);
        System.out.println("h_1_3 = " + h_1_3);
        System.out.println("-----------------------");
    }

    // gystota reshetki b/t = bt
    public double bt_1;
    public double A_1;
    public double B_1;
    public double C_1;
    public double J_1;

    public void gystota_pervoy_reshetki() {
        this.A_1 = otn_L_k_u_1 / otn_c_1_alfa_1;
        this.B_1 = Rou_st_1 / otn_c_1_alfa_1;
        this.C_1 = 0.7 - 0.27 * B_1 + 0.16 * B_1 * B_1;
        this.J_1 = A_1 / C_1;
        this.bt_1 = 0.225 + 0.275 * J_1 + 0.5 * J_1 * J_1;
        System.out.println("\n-----------------------");
        System.out.println("bt_1 = " + bt_1);
        System.out.println("A_1 = " + A_1);
        System.out.println("B_1 = " + B_1);
        System.out.println("C_1 = " + C_1);
        System.out.println("J_1 = " + J_1);
        System.out.println("-----------------------");
    }

    // gystota reshetki b/t = bt
    public double bt_2;
    public double A_2;
    public double B_2;
    public double C_2;
    public double J_2;

    public void gystota_vtorou_reshetki() {
        this.A_2 = otn_L_k_u_2 / otn_c_1_alfa_2;
        this.B_2 = Rou_st_2 / otn_c_1_alfa_2;
        this.C_2 = 0.7 - 0.27 * B_2 + 0.16 * B_2 * B_2;
        this.J_2 = A_2 / C_2;
        this.bt_2 = 0.225 + 0.275 * J_2 + 0.5 * J_2 * J_2;
        System.out.println("\n-----------------------");
        System.out.println("bt_2 = " + bt_2);
        System.out.println("A_2 = " + A_2);
        System.out.println("B_2 = " + B_2);
        System.out.println("C_2 = " + C_2);
        System.out.println("J_2 = " + J_2);
        System.out.println("-----------------------");
    }

    // gystota reshetki b/t = bt
    public double bt_3;
    public double A_3;
    public double B_3;
    public double C_3;
    public double J_3;

    public void gystota_tretiy_reshetki() {
        this.A_3 = otn_L_k_u_3 / otn_c_1_alfa_3;
        this.B_3 = Rou_st_3 / otn_c_1_alfa_3;
        this.C_3 = 0.7 - 0.27 * B_3 + 0.16 * B_3 * B_3;
        this.J_3 = A_3 / C_3;
        this.bt_3 = 0.225 + 0.275 * J_3 + 0.5 * J_3 * J_3;
        System.out.println("\n-----------------------");
        System.out.println("bt_3 = " + bt_3);
        System.out.println("A_3 = " + A_3);
        System.out.println("B_3 = " + B_3);
        System.out.println("C_3 = " + C_3);
        System.out.println("J_3 = " + J_3);
        System.out.println("-----------------------");
    }

    public double s_vna;
    public double s_rk_1;
    public double s_rk_2;
    public double s_rk_3;
    // эти параметры по рекомендации h_rk / s_rk
    public double h_rk_vna_del_s_rk_vna = 3.5;
    public double h_rk_1_del_s_rk_1 = 4.5;
    public double h_rk_2_del_s_rk_2 = 3.5;
    public double h_rk_3_del_s_rk_3 = 2.5;

    public void shirina_rab_lopatok_na_stypeniax() {
        this.s_vna = h_1_1 / h_rk_vna_del_s_rk_vna;
        this.s_rk_1 = h_1_1 / h_rk_1_del_s_rk_1;
        this.s_rk_2 = h_1_2 / h_rk_2_del_s_rk_2;
        this.s_rk_3 = h_1_3 / h_rk_3_del_s_rk_3;
        System.out.println("\n-----------------------");
        System.out.println("s_vna = " + s_vna);
        System.out.println("s_rk_1 = " + s_rk_1);
        System.out.println("s_rk_2 = " + s_rk_2);
        System.out.println("s_rk_3 = " + s_rk_3);
        System.out.println("-----------------------");
    }

    public double s_na_1;
    public double s_na_2;
    public double s_na_3;

    public void shirina_naprav_apparatov() {
        this.s_na_1 = 0.95 * s_rk_1;
        this.s_na_2 = 0.95 * s_rk_2;
        this.s_na_3 = 0.95 * s_rk_3;
        System.out.println("\n-----------------------");
        System.out.println("s_na_1 = " + s_na_1);
        System.out.println("s_na_2 = " + s_na_2);
        System.out.println("s_na_3 = " + s_na_3);
        System.out.println("-----------------------");
    }

    public double delta_a_1;
    public double delta_a_2;
    public double delta_a_3;

    public void osevie_zazori_mezhdy_venz_rab_koles_i_na() {
        this.delta_a_1 = 0.27 * s_rk_1;
        this.delta_a_2 = 0.27 * s_rk_2;
        this.delta_a_3 = 0.27 * s_rk_3;
        System.out.println("\n-----------------------");
        System.out.println("delta_a_1 = " + delta_a_1);
        System.out.println("delta_a_2 = " + delta_a_2);
        System.out.println("delta_a_3 = " + delta_a_3);
        System.out.println("-----------------------");
    }

    //poschitatb
    public double shirina_protoka;

    public void shirina_protochnoy_chasti() {
        System.out.println("shirina_protoka не посчитано");
    }

    public void calculate() {
        calculate_P_vx$();
        calculate_D_v_1();
        calculate_D_sr_1();
        calculate_D_vt_1();
        calculate_otn_F();
        calculate_F_v_vblx();
        calculate_F_v_vblx();
        calculate_otn_d_vblx();
        calculate_D_vblx_sr();
        calculate_D_vblx_vt();
        calculate_sum_otn_L_v();
        PODBOR_STYPENEY();
        calculate_L_v();
        calculate_sum_L_v();
        calculate_u_sr_1();
        koef_nagryzki_stypenei();
        calculate_N_v();
        calculate_F_vblx_1();
        calculate_F_vblx_2();
        calculate_D_razd();
        calculate_izontrop_potok();
        calculate_pi_v();
        polnoe_davlenie_na_vxode_v_stypeni();
        kriticheskie_skorosti_na_vxode();
        skorost_na_vxode_v_ventilyator();
        osevaya_skorost_kompressora();
        skorosti_na_stypenyax();
        koef_rasxoda_na_srednem_diametre_kolesa();
        koef_napora_na_kolese();
        ygol_cxods_vozdyxa_v_koleso_na_srednem_radiyse();
        prived_skorost_na_vxode();
        chislo_maxa_na_vxode_po_otn_skorosti();
        koef_proizvoditelnosti_pervoi_stypeni();
        okryzhnaya_sostavl_absolut_skorosti();
        absolut_i_prived_skorosti_na_vxode_v_rab_kolesa();
        ygol_vxoda_v_stypen_po_absolut_skorosti();
        calculate_ploshad_proxodnogo_sechenia();
        naryzniy_diametr_rab_koles();
        vnytrenniy_diametr_rab_koles();
        vblsota_rabochey_lopatki();
        gystota_pervoy_reshetki();
        gystota_vtorou_reshetki();
        gystota_tretiy_reshetki();
        shirina_rab_lopatok_na_stypeniax();
        shirina_naprav_apparatov();
        osevie_zazori_mezhdy_venz_rab_koles_i_na();
        shirina_protochnoy_chasti();
        var detalbniy = new P2_Детальный_расчет_вентилятора();
        detalbniy.calculate();
        P2_Детальный_расчет_вентилятора.Р3_Предварительный_расчет_турбины_низкого_давления predvarTND = detalbniy.new Р3_Предварительный_расчет_турбины_низкого_давления();
        predvarTND.calculate();
        var detalniyTND = predvarTND.new P4_Детальный_расчет_турбины_низкого_давления();
        detalniyTND.calculate();
        var paramPotoka = detalniyTND.new Р5_Расчет_параметров_потока_по_радиусу_ступени_турбины();
        paramPotoka.calculate();
    }

    class P2_Детальный_расчет_вентилятора {

        double betta_1_1;
        double betta_1_2;
        double betta_1_3;

        public void opredelenie_ygla_potoka_na_vxode_v_rk_v_otn() {
            double tg_betta_1_1 = otn_c_1_alfa_1 / (1 - otn_c_1_alfa_1 * ctg(alfa_1_1_f));
            double tg_betta_1_2 = otn_c_1_alfa_2 / (1 - otn_c_1_alfa_2 * ctg(alfa_1_2_f));
            double tg_betta_1_3 = otn_c_1_alfa_3 / (1 - otn_c_1_alfa_3 * ctg(alfa_1_3_f));
            System.out.println("\n-----------------------");
            System.out.println("ctg(alfa_1_1_f) = " + ctg(alfa_1_1_f));
            System.out.println("ctg(alfa_1_2_f) = " + ctg(alfa_1_2_f));
            System.out.println("ctg(alfa_1_3_f) = " + ctg(alfa_1_3_f));
            System.out.println("tg_betta_1_1 = " + tg_betta_1_1);
            System.out.println("tg_betta_1_2 = " + tg_betta_1_2);
            System.out.println("tg_betta_1_3 = " + tg_betta_1_3);
            this.betta_1_1 = fromRadToGrad(Math.atan(tg_betta_1_1));
            this.betta_1_2 = fromRadToGrad(Math.atan(tg_betta_1_2));
            this.betta_1_3 = fromRadToGrad(Math.atan(tg_betta_1_3));
            System.out.println("betta_1_1 = " + betta_1_1);
            System.out.println("betta_1_2 = " + betta_1_2);
            System.out.println("betta_1_3 = " + betta_1_3);
            System.out.println("-----------------------");
        }

        double W_1_1;
        double W_1_2;
        double W_1_3;

        public void otn_skorosti_na_vxode_v_RK() {
            this.W_1_1 = Math.sqrt(c_1_a_1 * c_1_a_1 + (u_sr_1 - c_1_u_1) * (u_sr_1 - c_1_u_1));
            this.W_1_2 = Math.sqrt(c_1_a_2 * c_1_a_2 + (u_sr_2 - c_1_u_2) * (u_sr_2 - c_1_u_2));
            this.W_1_3 = Math.sqrt(c_1_a_3 * c_1_a_3 + (u_sr_3 - c_1_u_3) * (u_sr_3 - c_1_u_3));
            System.out.println("\n-----------------------");
            System.out.println("W_1_1 = " + W_1_1);
            System.out.println("W_1_2 = " + W_1_2);
            System.out.println("W_1_3 = " + W_1_3);
            System.out.println("-----------------------");
        }

        double a_1_1;
        double a_1_2;
        double a_1_3;
        double a_1_vblx;
        // смотрим на lambda_1_1_f
        double tau_lambda_1_1 = 0.9229;
        double tau_lambda_1_2 = 0.9317;
        double tau_lambda_1_3 = 0.9420;
        double tau_lambda_v_vblx = 0.9496;

        public void skorosti_na_vxode() {
            this.a_1_1 = alfa_1_kr_1 * Math.sqrt((k_vozdyxa + 1) / 2 * tau_lambda_1_1);
            this.a_1_2 = alfa_1_kr_2 * Math.sqrt((k_vozdyxa + 1) / 2 * tau_lambda_1_2);
            this.a_1_3 = alfa_1_kr_3 * Math.sqrt((k_vozdyxa + 1) / 2 * tau_lambda_1_3);
            this.a_1_vblx = alfa_1_kr_vblx * Math.sqrt((k_vozdyxa + 1) / 2 * tau_lambda_v_vblx);
            System.out.println("\n-----------------------");
            System.out.println("a_1_1 = " + a_1_1);
            System.out.println("a_1_2 = " + a_1_2);
            System.out.println("a_1_3 = " + a_1_3);
            System.out.println("a_1_vblx = " + a_1_vblx);
            System.out.println("-----------------------");
        }

        double M_w_1_1;
        double M_w_1_2;
        double M_w_1_3;

        public void chislo_Maxa_v_otn_dv_na_vxode_v_stypeni() {
            this.M_w_1_1 = W_1_1 / a_1_1;
            this.M_w_1_2 = W_1_2 / a_1_2;
            this.M_w_1_3 = W_1_3 / a_1_3;
            System.out.println("\n-----------------------");
            System.out.println("M_w_1_1 = " + M_w_1_1);
            System.out.println("M_w_1_2 = " + M_w_1_2);
            System.out.println("M_w_1_3 = " + M_w_1_3);
            System.out.println("-----------------------");
        }

        double h_ydl_1 = 3.71;
        double h_ydl_2 = 3.48;
        double h_ydl_3 = 3.26;

        int z_rk_1;
        int z_rk_2;
        int z_rk_3;

        public void chislo_lopatok_rab_koles() {
            System.out.println("\n-----------------------");
            System.out.println("h_ydl_1 [ДО] = " + h_ydl_1);
            z_rk_1 = (int) (h_ydl_1 * bt_1 * Math.PI * D_sr_1 / h_1_1);
            System.out.println("z_rk_1 = " + z_rk_1);
            h_ydl_1 = (z_rk_1 * h_1_1) / (bt_1 * Math.PI * D_sr_1);
            System.out.println("h_ydl_1 [ПОСЛЕ] = " + h_ydl_1);
            System.out.println("-----------------------");

            System.out.println("\n-----------------------");
            System.out.println("h_ydl_2 [ДО] = " + h_ydl_2);
            z_rk_2 = (int) (h_ydl_2 * bt_2 * Math.PI * D_sr_2 / h_1_2);
            System.out.println("z_rk_2 = " + z_rk_2);
            h_ydl_2 = (z_rk_2 * h_1_2) / (bt_2 * Math.PI * D_sr_2);
            System.out.println("h_ydl_2 [ПОСЛЕ] = " + h_ydl_2);
            System.out.println("-----------------------");

            System.out.println("\n-----------------------");
            System.out.println("h_ydl_3 [ДО] = " + h_ydl_3);
            z_rk_3 = (int) (h_ydl_3 * bt_3 * Math.PI * D_sr_3 / h_1_3);
            System.out.println("z_rk_3 = " + z_rk_3);
            h_ydl_3 = (z_rk_3 * h_1_3) / (bt_3 * Math.PI * D_sr_3);
            System.out.println("h_ydl_3 [ПОСЛЕ] = " + h_ydl_3);
            System.out.println("-----------------------");
        }

        //хорды
        double b_rk_1;
        double b_rk_2;
        double b_rk_3;

        public void dlina_xord_lopatok_ventilyatora() {
            this.b_rk_1 = h_1_1 / h_ydl_1;
            this.b_rk_2 = h_1_2 / h_ydl_2;
            this.b_rk_3 = h_1_3 / h_ydl_3;
            System.out.println("\n-----------------------");
            System.out.println("b_rk_1 = " + b_rk_1);
            System.out.println("b_rk_2 = " + b_rk_2);
            System.out.println("b_rk_3 = " + b_rk_3);
            System.out.println("-----------------------");
        }

        double c_2_u_1;
        double c_2_u_2;
        double c_2_u_3;

        public void okr_sostavl_absolutnoy_skorosti_na_vblxpde_iz_rab_kolesa() {
            this.c_2_u_1 = u_sr_1 * ((1 - Rou_st_1) + (otn_L_k_u_1 / 2));
            this.c_2_u_2 = u_sr_2 * ((1 - Rou_st_2) + (otn_L_k_u_2 / 2));
            this.c_2_u_3 = u_sr_3 * ((1 - Rou_st_3) + (otn_L_k_u_3 / 2));
            System.out.println("\n-----------------------");
            System.out.println("c_2_u_1 = " + c_2_u_1);
            System.out.println("c_2_u_2 = " + c_2_u_2);
            System.out.println("c_2_u_3 = " + c_2_u_3);
            System.out.println("-----------------------");
        }

        double c_2_a_1;
        double c_2_a_2;
        double c_2_a_3;

        public void osevie_skorosti_na_vblxode_iz_rab_kolesa_stypeniy() {
            this.c_2_a_1 = (c_1_a_1 + c_1_a_2) / 2;
            this.c_2_a_2 = (c_1_a_2 + c_1_a_3) / 2;
            this.c_2_a_3 = (c_1_a_3 + c_1_a_vblx) / 2;
            System.out.println("\n-----------------------");
            System.out.println("c_2_a_1 = " + c_2_a_1);
            System.out.println("c_2_a_2 = " + c_2_a_2);
            System.out.println("c_2_a_3 = " + c_2_a_3);
            System.out.println("-----------------------");
        }

        double c_2_1;
        double c_2_2;
        double c_2_3;
        double lambda_2_1;
        double lambda_2_2;
        double lambda_2_3;

        public void absol_i_prived_skorost_na_vblxode_iz_rab_koles_stypeney_ventilyatora() {
            this.c_2_1 = Math.sqrt(c_2_a_1 * c_2_a_1 + c_2_u_1 * c_2_u_1);
            this.lambda_2_1 = c_2_1 / a_1_1;
            this.c_2_2 = Math.sqrt(c_2_a_2 * c_2_a_2 + c_2_u_2 * c_2_u_2);
            this.lambda_2_2 = c_2_2 / a_1_2;
            this.c_2_3 = Math.sqrt(c_2_a_3 * c_2_a_3 + c_2_u_3 * c_2_u_3);
            this.lambda_2_3 = c_2_3 / a_1_3;
            System.out.println("\n-----------------------");
            System.out.println("c_2_1 = " + c_2_1);
            System.out.println("lambda_2_1 = " + lambda_2_1);
            System.out.println("c_2_2 = " + c_2_2);
            System.out.println("lambda_2_2 = " + lambda_2_2);
            System.out.println("c_2_3 = " + c_2_3);
            System.out.println("lambda_2_3 = " + lambda_2_3);
            System.out.println("-----------------------");
        }

        double a_2_1;
        double a_2_2;
        double a_2_3;
        // смотрим на lambda_1_1_f
        double tau_lambda_2_1 = 0.8739;
        double tau_lambda_2_2 = 0.8824;
        double tau_lambda_2_3 = 0.8986;

        public void mestniy_skorosti_na_vblxode_iz_rabochix_koles() {
            this.a_2_1 = alfa_1_kr_2 * Math.sqrt((k_vozdyxa + 1) / 2 * tau_lambda_2_1);
            this.a_2_2 = alfa_1_kr_3 * Math.sqrt((k_vozdyxa + 1) / 2 * tau_lambda_2_2);
            this.a_2_3 = alfa_1_kr_vblx * Math.sqrt((k_vozdyxa + 1) / 2 * tau_lambda_2_3);
            System.out.println("\n-----------------------");
            System.out.println("a_2_1 = " + a_2_1);
            System.out.println("a_2_2 = " + a_2_2);
            System.out.println("a_2_3 = " + a_2_3);
            System.out.println("-----------------------");
        }

        double M_c_2_1;
        double M_c_2_2;
        double M_c_2_3;

        public void chislo_Maxa_po_absolut_skorosti_na_vxode_v_napravl_app_stypeney_ventilyatora() {
            this.M_c_2_1 = c_2_1 / a_2_1;
            this.M_c_2_2 = c_2_2 / a_2_2;
            this.M_c_2_3 = c_2_3 / a_2_3;
            System.out.println("\n-----------------------");
            System.out.println("M_c_2_1 = " + M_c_2_1);
            System.out.println("M_c_2_2 = " + M_c_2_2);
            System.out.println("M_c_2_3 = " + M_c_2_3);
            System.out.println("-----------------------");
        }

        double sin_alfa_2_1;
        double sin_alfa_2_2;
        double sin_alfa_2_3;
        double alfa_2_1;
        double alfa_2_2;
        double alfa_2_3;

        public void opred_yglov_vblxoda_iz_rab_koles_v_absolut_dvizhenii() {
            this.sin_alfa_2_1 = c_2_a_1 / c_2_1;
            this.sin_alfa_2_2 = c_2_a_2 / c_2_2;
            this.sin_alfa_2_3 = c_2_a_3 / c_2_3;
            this.alfa_2_1 = fromRadToGrad(Math.asin(sin_alfa_2_1));
            this.alfa_2_2 = fromRadToGrad(Math.asin(sin_alfa_2_2));
            this.alfa_2_3 = fromRadToGrad(Math.asin(sin_alfa_2_3));
            System.out.println("\n-----------------------");
            System.out.println("sin_alfa_2_1 = " + sin_alfa_2_1);
            System.out.println("sin_alfa_2_2 = " + sin_alfa_2_2);
            System.out.println("sin_alfa_2_3 = " + sin_alfa_2_3);
            System.out.println("alfa_2_1 = " + alfa_2_1);
            System.out.println("alfa_2_2 = " + alfa_2_2);
            System.out.println("alfa_2_3 = " + alfa_2_3);
            System.out.println("-----------------------");
        }

        double P_2_1$;
        double P_2_2$;
        double P_2_3$;

        public void poln_davl_na_vblxode_iz_rabochix_koles_stypenei_ventilyatora() {
            this.P_2_1$ = P_1_1$ * Math.pow((1 + ((L_v_1 * Nu_rk) / (k_vozdyxa / (k_vozdyxa - 1) * R_vozdyxa * T_vx$))), k_vozdyxa / (k_vozdyxa - 1));
            this.P_2_2$ = P_1_2$ * Math.pow((1 + ((L_v_2 * Nu_rk) / (k_vozdyxa / (k_vozdyxa - 1) * R_vozdyxa * T_1_2$))), k_vozdyxa / (k_vozdyxa - 1));
            this.P_2_3$ = P_1_3$ * Math.pow((1 + ((L_v_3 * Nu_rk) / (k_vozdyxa / (k_vozdyxa - 1) * R_vozdyxa * T_1_3$))), k_vozdyxa / (k_vozdyxa - 1));
            System.out.println("\n-----------------------");
            System.out.println("P_2_1$ = " + P_2_1$);
            System.out.println("P_2_2$ = " + P_2_2$);
            System.out.println("P_2_3$ = " + P_2_3$);
            System.out.println("-----------------------");
        }

        double sigma_CA_1;
        double sigma_CA_2;
        double sigma_CA_3;

        public void koef_vostan_polnogo_davleniya() {
            this.sigma_CA_1 = P_1_2$ / P_2_1$;
            this.sigma_CA_2 = P_1_3$ / P_2_2$;
            this.sigma_CA_3 = P_v_vblx$ / P_2_3$;
            System.out.println("\n-----------------------");
            System.out.println("sigma_CA_1 = " + sigma_CA_1);
            System.out.println("sigma_CA_2 = " + sigma_CA_2);
            System.out.println("sigma_CA_3 = " + sigma_CA_3);
            System.out.println("-----------------------");
        }

        double F_2_1;
        double F_2_2;
        double F_2_3;

        double q_lambda_2_1 = 0.9796;
        double q_lambda_2_2 = 0.9691;
        double q_lambda_2_3 = 0.9418;

        public void calculate_ploshad_proxodnogo_sechenia() {
            this.F_2_1 = (G_v_sum * Math.sqrt(T_2_1$)) /
                    (P_2_1$ * q_lambda_2_1 * sin_alfa_2_1 * s_v * Kg);
            this.F_2_2 = (G_v_sum * Math.sqrt(T_2_2$)) /
                    (P_2_2$ * q_lambda_2_2 * sin_alfa_2_2 * s_v * Kg);
            this.F_2_3 = (G_v_sum * Math.sqrt(T_3_3$)) /
                    (P_2_3$ * q_lambda_2_3 * sin_alfa_2_3 * s_v * Kg);
            System.out.println("\n-----------------------");
            System.out.println("q_lambda_2_1 = " + q_lambda_2_1);
            System.out.println("q_lambda_2_2 = " + q_lambda_2_2);
            System.out.println("q_lambda_2_3 = " + q_lambda_2_3);
            System.out.println("F_2_1 = " + F_2_1);
            System.out.println("F_2_2 = " + F_2_2);
            System.out.println("F_2_3 = " + F_2_3);
            System.out.println("-----------------------");
        }

        double otn_d_vt_2_1;
        double otn_d_vt_2_2;
        double otn_d_vt_2_3;

        public void otn_diametr_vtylki_za_rab_kolesami_stypeney() {
            this.otn_d_vt_2_1 = Math.sqrt(1 - (4 * F_2_1) / (Math.PI * D_v * D_v));
            this.otn_d_vt_2_2 = Math.sqrt(1 - (4 * F_2_2) / (Math.PI * D_v * D_v));
            this.otn_d_vt_2_3 = Math.sqrt(1 - (4 * F_2_3) / (Math.PI * D_v * D_v));
            System.out.println("\n-----------------------");
            System.out.println("otn_d_vt_2_1 = " + otn_d_vt_2_1);
            System.out.println("otn_d_vt_2_2 = " + otn_d_vt_2_2);
            System.out.println("otn_d_vt_2_3 = " + otn_d_vt_2_3);
            System.out.println("-----------------------");
        }

        double D_vt_2_1;
        double D_vt_2_2;
        double D_vt_2_3;

        public void diametr_vtylki_za_rab_kolesami_stypeney_ventilyatora() {
            this.D_vt_2_1 = D_v * otn_d_vt_2_1;
            this.D_vt_2_2 = D_v * otn_d_vt_2_2;
            this.D_vt_2_3 = D_v * otn_d_vt_2_3;
            System.out.println("\n-----------------------");
            System.out.println("D_2_vt_1 = " + D_vt_2_1);
            System.out.println("D_2_vt_2 = " + D_vt_2_2);
            System.out.println("D_2_vt_3 = " + D_vt_2_3);
            System.out.println("-----------------------");
        }

        double h_2_1;
        double h_2_2;
        double h_2_3;

        public void visoti_lopatok_na_vblxode_iz_rab_koles_stypenei_ventilyatora() {
            this.h_2_1 = (D_v - D_vt_1_1) / 2;
            this.h_2_2 = (D_v - D_vt_1_2) / 2;
            this.h_2_3 = (D_v - D_vt_1_3) / 2;
            System.out.println("\n-----------------------");
            System.out.println("h_2_1 = " + h_2_1);
            System.out.println("h_2_2 = " + h_2_2);
            System.out.println("h_2_3 = " + h_2_3);
            System.out.println("-----------------------");
        }

        double W_2_1;
        double W_2_2;
        double W_2_3;

        public void otn_skorosti_na_vblxode_iz_RK() {
            this.W_2_1 = Math.sqrt(c_2_a_1 * c_2_a_1 + (u_sr_1 - c_2_u_1) * (u_sr_1 - c_2_u_1));
            this.W_2_2 = Math.sqrt(c_2_a_2 * c_2_a_2 + (u_sr_2 - c_2_u_2) * (u_sr_2 - c_2_u_2));
            this.W_2_3 = Math.sqrt(c_2_a_3 * c_2_a_3 + (u_sr_3 - c_2_u_3) * (u_sr_3 - c_2_u_3));
            System.out.println("\n-----------------------");
            System.out.println("W_2_1 = " + W_2_1);
            System.out.println("W_2_2 = " + W_2_2);
            System.out.println("W_2_3 = " + W_2_3);
            System.out.println("-----------------------");
        }

        double sin_betta_2_1;
        double sin_betta_2_2;
        double sin_betta_2_3;
        double betta_2_1;
        double betta_2_2;
        double betta_2_3;

        public void ygli_iz_rabochix_koles_v_otn_dvizhenii() {
            this.sin_betta_2_1 = c_2_a_1 / W_2_1;
            this.sin_betta_2_2 = c_2_a_2 / W_2_2;
            this.sin_betta_2_3 = c_2_a_3 / W_2_3;
            this.betta_2_1 = fromRadToGrad(Math.asin(sin_betta_2_1));
            this.betta_2_2 = fromRadToGrad(Math.asin(sin_betta_2_2));
            this.betta_2_3 = fromRadToGrad(Math.asin(sin_betta_2_3));
            System.out.println("\n-----------------------");
            System.out.println("sin_betta_2_1 = " + sin_betta_2_1);
            System.out.println("sin_betta_2_2 = " + sin_betta_2_2);
            System.out.println("sin_betta_2_3 = " + sin_betta_2_3);
            System.out.println("betta_2_1 = " + betta_2_1);
            System.out.println("betta_2_2 = " + betta_2_2);
            System.out.println("betta_2_3 = " + betta_2_3);
            System.out.println("-----------------------");
        }

        double delta_betta_1;
        double delta_betta_2;
        double delta_betta_3;

        public void ygli_povorota_v_rab_kolesax_styp_ventilyatora() {
            this.delta_betta_1 = betta_2_1 - betta_1_1;
            this.delta_betta_2 = betta_2_2 - betta_1_2;
            this.delta_betta_3 = betta_2_3 - betta_1_3;
            System.out.println("\n-----------------------");
            System.out.println("delta_betta_1 = " + delta_betta_1);
            System.out.println("delta_betta_2 = " + delta_betta_2);
            System.out.println("delta_betta_3 = " + delta_betta_3);
            System.out.println("-----------------------");
        }

        double alfa_3_1;
        double alfa_3_2;
        double alfa_3_vblx;

        public void ygli_vblxoda_potoka_iz_napravl_app() {
            this.alfa_3_1 = alfa_1_2_f;
            this.alfa_3_2 = alfa_1_3_f;
            this.alfa_3_vblx = 90;

            System.out.println("\n-----------------------");
            System.out.println("alfa_3_1 = " + alfa_3_1);
            System.out.println("alfa_3_2 = " + alfa_3_2);
            System.out.println("alfa_3_vblx = " + alfa_3_vblx);
            System.out.println("-----------------------");
        }

        double delta_alfa_1;
        double delta_alfa_2;
        double delta_alfa_3;

        public void ygli_povorota_potoka_v_napravl_app_stypeney_ventilyatora() {
            this.delta_alfa_1 = alfa_3_1 - alfa_2_1;
            this.delta_alfa_2 = alfa_3_2 - alfa_2_2;
            this.delta_alfa_3 = alfa_3_vblx - alfa_2_3;
            System.out.println("\n-----------------------");
            System.out.println("delta_alfa_1 = " + delta_alfa_1);
            System.out.println("delta_alfa_2 = " + delta_alfa_2);
            System.out.println("delta_alfa_3 = " + delta_alfa_3);
            System.out.println("-----------------------");
        }

        double otn_alfa_3_1;
        double otn_alfa_3_2;
        double otn_alfa_3_3;
        double otn_delta_alfa_1_1;
        double otn_delta_alfa_1_2;
        double otn_delta_alfa_1_3;
        double otn_delta_alfa_1;
        double otn_delta_alfa_2;
        double otn_delta_alfa_3;
        double E_1;
        double E_2;
        double E_3;

        public void nominalniy_egol_povorota_v_absolut_dvizh_na_stypen_vent() {
            System.out.println("\n-----------------------");
            //номинальный угол потока на выходе из ступени вентилятора
            this.otn_alfa_3_1 = alfa_3_1 / 100;
            this.otn_alfa_3_2 = alfa_3_2 / 100;
            this.otn_alfa_3_3 = alfa_3_vblx / 100;
            System.out.println("otn_alfa_3_1 = " + otn_alfa_3_1);
            System.out.println("otn_alfa_3_2 = " + otn_alfa_3_2);
            System.out.println("otn_alfa_3_3 = " + otn_alfa_3_3);

            //номинальные углы поворота потока за ступенями вентилятора, при густоте решетки равной единице
            this.otn_delta_alfa_1_1 = 0.037 + 0.1 * otn_alfa_3_1 + 0.262 * otn_alfa_3_1 * otn_alfa_3_1;
            this.otn_delta_alfa_1_2 = 0.037 + 0.1 * otn_alfa_3_2 + 0.262 * otn_alfa_3_2 * otn_alfa_3_2;
            this.otn_delta_alfa_1_3 = 0.037 + 0.1 * otn_alfa_3_3 + 0.262 * otn_alfa_3_3 * otn_alfa_3_3;
            System.out.println("otn_delta_alfa_1_1 = " + otn_delta_alfa_1_1);
            System.out.println("otn_delta_alfa_1_2 = " + otn_delta_alfa_1_2);
            System.out.println("otn_delta_alfa_1_3 = " + otn_delta_alfa_1_3);

            this.otn_delta_alfa_1 = delta_alfa_1 / 100;
            this.otn_delta_alfa_2 = delta_alfa_2 / 100;
            this.otn_delta_alfa_3 = delta_alfa_3 / 100;
            System.out.println("otn_delta_alfa_1 = " + otn_delta_alfa_1);
            System.out.println("otn_delta_alfa_2 = " + otn_delta_alfa_2);
            System.out.println("otn_delta_alfa_3 = " + otn_delta_alfa_3);

            //определяется параметр E
            this.E_1 = otn_delta_alfa_1 / otn_delta_alfa_1_1;
            this.E_2 = otn_delta_alfa_2 / otn_delta_alfa_1_2;
            this.E_3 = otn_delta_alfa_3 / otn_delta_alfa_1_3;
            System.out.println("E_1 = " + E_1);
            System.out.println("E_2 = " + E_2);
            System.out.println("E_3 = " + E_3);
            System.out.println("-----------------------");
        }

        double bt_CA_1;
        double bt_CA_2;
        double bt_CA_3;

        //Определение густоты решетки спрямляющих аппаратов на ступенях вентилятора
        public void opred_gystoti_resh_v_spryamlyash_apparatax() {
            this.bt_CA_1 = 0.231 - 0.135 * E_1 + 0.909 * E_1 * E_1;
            this.bt_CA_2 = 10 * (0.981 - 1.788 * E_2 + 0.912 * E_2 * E_2);
            this.bt_CA_3 = 10 * (0.981 - 1.788 * E_3 + 0.912 * E_3 * E_3);
            System.out.println("\n-----------------------");
            System.out.println("bt_CA_1 = " + bt_CA_1);
            System.out.println("bt_CA_2 = " + bt_CA_2);
            System.out.println("bt_CA_3 = " + bt_CA_3);
            System.out.println("-----------------------");
        }

        int z_CA_1;
        int z_CA_2;
        int z_CA_3;
        double h_ydl_CA_1;
        double h_ydl_CA_2;
        double h_ydl_CA_3;

        //Алгоритм расчета такой же, как и для определения числа лопаток рабочих колес. Удлинения лопаток будут взяты такие же, как и для рабочих колес. Тогда
        public void opredelenie_chisla_lopatok_spryamlyash_apparatov() {
            this.h_ydl_CA_1 = h_ydl_1;
            this.h_ydl_CA_2 = h_ydl_2;
            this.h_ydl_CA_3 = h_ydl_3;
            System.out.println("\n-----------------------");
            System.out.println("h_ydl_CA_1 [ДО] = " + h_ydl_CA_1);
            this.z_CA_1 = (int) ((h_ydl_CA_1 * bt_CA_1 * Math.PI * D_sr_1) / h_2_1);
            System.out.println("z_CA_1 = " + z_CA_1);
            h_ydl_CA_1 = (z_CA_1 * h_2_1) / (bt_CA_1 * Math.PI * D_sr_1);
            System.out.println("h_ydl_CA_1 [ПОСЛЕ] = " + h_ydl_CA_1);
            System.out.println("-----------------------");

            System.out.println("\n-----------------------");
            System.out.println("h_ydl_CA_2 [ДО] = " + h_ydl_CA_2);
            this.z_CA_2 = (int) ((h_ydl_CA_2 * bt_CA_2 * Math.PI * D_sr_2) / h_2_2);
            System.out.println("z_CA_2 = " + z_CA_2);
            h_ydl_CA_2 = (z_CA_2 * h_2_2) / (bt_CA_2 * Math.PI * D_sr_2);
            System.out.println("h_ydl_2 [ПОСЛЕ] = " + h_ydl_CA_2);
            System.out.println("-----------------------");

            System.out.println("\n-----------------------");
            System.out.println("h_ydl_CA_3 [ДО] = " + h_ydl_CA_3);
            this.z_CA_3 = (int) ((h_ydl_CA_3 * bt_CA_3 * Math.PI * D_sr_3) / h_2_3);
            System.out.println("z_CA_3 = " + z_CA_3);
            h_ydl_CA_3 = (z_CA_3 * h_2_3) / (bt_CA_3 * Math.PI * D_sr_3);
            System.out.println("h_ydl_3 [ПОСЛЕ] = " + h_ydl_CA_3);
            System.out.println("-----------------------");
        }

        double b_CA_1;
        double b_CA_2;
        double b_CA_3;

        public void opred_hord_streml_apparatov() {
            this.b_CA_1 = h_2_1 / h_ydl_CA_1;
            this.b_CA_2 = h_2_2 / h_ydl_CA_2;
            this.b_CA_3 = h_2_3 / h_ydl_CA_3;
            System.out.println("\n-----------------------");
            System.out.println("b_CA_1 = " + b_CA_1);
            System.out.println("b_CA_2 = " + b_CA_2);
            System.out.println("b_CA_3 = " + b_CA_3);
            System.out.println("-----------------------");
        }

        public void calculate() {
            opredelenie_ygla_potoka_na_vxode_v_rk_v_otn();
            otn_skorosti_na_vxode_v_RK();
            skorosti_na_vxode();
            chislo_Maxa_v_otn_dv_na_vxode_v_stypeni();
            chislo_lopatok_rab_koles();
            dlina_xord_lopatok_ventilyatora();
            okr_sostavl_absolutnoy_skorosti_na_vblxpde_iz_rab_kolesa();
            osevie_skorosti_na_vblxode_iz_rab_kolesa_stypeniy();
            absol_i_prived_skorost_na_vblxode_iz_rab_koles_stypeney_ventilyatora();
            mestniy_skorosti_na_vblxode_iz_rabochix_koles();
            chislo_Maxa_po_absolut_skorosti_na_vxode_v_napravl_app_stypeney_ventilyatora();
            opred_yglov_vblxoda_iz_rab_koles_v_absolut_dvizhenii();
            poln_davl_na_vblxode_iz_rabochix_koles_stypenei_ventilyatora();
            koef_vostan_polnogo_davleniya();
            calculate_ploshad_proxodnogo_sechenia();
            otn_diametr_vtylki_za_rab_kolesami_stypeney();
            diametr_vtylki_za_rab_kolesami_stypeney_ventilyatora();
            visoti_lopatok_na_vblxode_iz_rab_koles_stypenei_ventilyatora();
            otn_skorosti_na_vblxode_iz_RK();
            ygli_iz_rabochix_koles_v_otn_dvizhenii();
            ygli_povorota_v_rab_kolesax_styp_ventilyatora();
            ygli_vblxoda_potoka_iz_napravl_app();
            ygli_povorota_potoka_v_napravl_app_stypeney_ventilyatora();
            nominalniy_egol_povorota_v_absolut_dvizh_na_stypen_vent();
            opred_gystoti_resh_v_spryamlyash_apparatax();
            opredelenie_chisla_lopatok_spryamlyash_apparatov();
            opred_hord_streml_apparatov();
        }

        class Р3_Предварительный_расчет_турбины_низкого_давления {
            //Так как расчитываемый двигатель, имеет не высокую степень двухконтурности, то данное отношение следует принять равным
            double D_tv_vblx_nar_del_D_razd = 0.93;
            double D_tv_vblx_nar;

            public void nar_diametr_tyrbinbl_na_vblxode() {
                this.D_tv_vblx_nar = D_tv_vblx_nar_del_D_razd * D_razd;
                System.out.println("\n-----------------------");
                System.out.println("D_tv_vblx_nar = " + D_tv_vblx_nar);
                System.out.println("-----------------------");
            }

            double L_tv;
            double k_t = 0.99;

            public void potrebnaya_vnytr_ydelnaya_rabota_tyrb_tyrboventilyatora() {
                this.L_tv = ((m + 1) * L_v) / k_t;
                System.out.println("\n-----------------------");
                System.out.println("L_tv = " + L_tv);
                System.out.println("-----------------------");
            }

            double T_tv$;
            double L_kvd;
            double N_kvd;
            double N_tvd;
            double L_tk;
            double T_vblx_tk$;

            public void temperatyra_gaza_pered_tyrbinoy() {
                this.L_kvd = 1005 * T_k_vx$ * (Math.pow(Pi_kvd$, (k_gaza - 1) / k_gaza) - 1) * (1 / Nu_kvd);
                this.N_kvd = L_kvd * G_v_1;
                this.N_tvd = N_kvd / Nu_tvd;
                this.L_tk = N_tvd / G_v_1;
                this.T_vblx_tk$ = T_g$ - (L_tk) / (k_gaza / (k_gaza - 1) * R_gaza);
                this.T_tv$ = T_vblx_tk$ - (L_tv) / (k_gaza / (k_gaza - 1) * R_gaza);

                System.out.println("\n-----------------------");
                System.out.println("L_kvd = " + L_kvd);
                System.out.println("N_kvd = " + N_kvd);
                System.out.println("N_tvd = " + N_tvd);
                System.out.println("L_tk = " + L_tk);
                System.out.println("T_vblx_tk$ = " + T_vblx_tk$);
                System.out.println("T_tv$ = " + T_tv$);
                System.out.println("-----------------------");
            }

            //Отношение полных давлений в турбине вентилятора
            double P_tk$_sigma_per_kan_del_P_tv$;

            public void calculate_poln_davleniy_v_tyrbine_ventilyatora() {
                double preDenominator = (k_gaza / (k_gaza - 1)) * R_gaza * T_vblx_tk$ * Nu_tv$;
                double denominator = Math.pow(1 - L_tv / preDenominator, k_gaza / (k_gaza - 1));
                this.P_tk$_sigma_per_kan_del_P_tv$ = 1 / denominator;
                System.out.println("\n-----------------------");
                System.out.println("P_tk$_sigma_per_kan_del_P_tv$ = " + P_tk$_sigma_per_kan_del_P_tv$);
                System.out.println("-----------------------");
            }

            double F_2_tv;
            double G_t;
            double sigma_per_kan = 0.99;
            double P_g$_del_P_tk$;
            double P_tk$;
            double P_tv$;
            double sin_alfa_2_tv;
            double q_lambda_2_tv = 0.6896;

            public void ploshad_kol_sechenia_na_vblxode_iz_tyrbinbl_ventilyatora() {
                this.G_t = 0.95 * G_v_1;
                double preDenominator = (k_gaza / (k_gaza - 1)) * R_gaza * T_g$ * Nu_tvd;
                double denominator = Math.pow(1 - L_tk / preDenominator, k_gaza / (k_gaza - 1));
                this.P_g$_del_P_tk$ = 1 / denominator;
                this.P_tk$ = (P_v_vblx$ * Pi_kvd$ * 0.97) / P_g$_del_P_tk$;
                this.P_tv$ = (P_tk$ * sigma_per_kan) / P_tk$_sigma_per_kan_del_P_tv$;
                this.sin_alfa_2_tv = Math.sin(fromGradToRad(alfa_2_t));
                this.F_2_tv = (G_t * Math.sqrt(T_tv$)) / (s_g * P_tv$ * q_lambda_2_tv * sin_alfa_2_tv);
                System.out.println("\n-----------------------");
                System.out.println("G_t = " + G_t);
                System.out.println("sigma_per_kan = " + sigma_per_kan);
                System.out.println("P_g$_del_P_tk$ = " + P_g$_del_P_tk$);
                System.out.println("P_tk$ = " + P_tk$);
                System.out.println("P_tv$ = " + P_tv$);
                System.out.println("sin_alfa_2_t = " + sin_alfa_2_tv);
                System.out.println("q_lambda_2_tv = " + q_lambda_2_tv);
                System.out.println("F_2_tv = " + F_2_tv);
                System.out.println("-----------------------");
            }

            double h_2_tv;

            public void vblsota_lopatok_tyrbinbl_vent() {
                double underSqrt = (D_tv_vblx_nar * D_tv_vblx_nar / 4) - (F_2_tv / Math.PI);
                h_2_tv = (D_tv_vblx_nar / 2) - Math.sqrt(underSqrt);
                System.out.println("\n-----------------------");
                System.out.println("h_2_tv = " + h_2_tv);
                System.out.println("-----------------------");
            }

            double D_tv_vblx_sr;

            public void sredniy_diameter_tyrbinbl_na_vblxode() {
                this.D_tv_vblx_sr = Math.sqrt((D_tv_vblx_nar * D_tv_vblx_nar) - (2 * F_2_tv) / Math.PI);
                System.out.println("\n-----------------------");
                System.out.println("D_tv_vblx_sr = " + D_tv_vblx_sr);
                System.out.println("-----------------------");
            }

            // D_sr = const m < 5

            double F_1_tv;
            double lambda_1_tv = 0.47;
            double q_lambda_1_tv = 0.6780;

            public void ploshad_protochnoy_chasti_na_vxode_v_tyrbiny_ventilyatora() {
                this.F_1_tv = (G_t * Math.sqrt(T_vblx_tk$)) / (s_g * P_tk$ * sigma_per_kan * q_lambda_1_tv);
                System.out.println("\n-----------------------");
                System.out.println("F_1_tv = " + F_1_tv);
                System.out.println("-----------------------");
            }

            double D_tv_sr;
            double D_1_tv_nar;
            double D_1_tv_vt;

            public void diametraln_razmeri_na_vxode_v_tyrbiny() {
                this.D_tv_sr = D_tv_vblx_sr;
                this.D_1_tv_nar = D_tv_sr + F_1_tv / (Math.PI * D_tv_sr);
                this.D_1_tv_vt = D_tv_sr - F_1_tv / (Math.PI * D_tv_sr);
                System.out.println("\n-----------------------");
                System.out.println("D_tv_sr = " + D_tv_sr);
                System.out.println("D_1_tv_nar = " + D_1_tv_nar);
                System.out.println("D_1_tv_vt = " + D_1_tv_vt);
                System.out.println("-----------------------");
            }

            // в первом приближении
            double K_tv = 0.45;
            // chislo stypeney ventilyatora
            double z_v = 3;
            int z_tv;

            public void chislo_stypeney_tyrbinbl() {
                this.z_tv = (int) ((K_tv * K_tv * z_v * (m + 1)) / ((D_tv_sr / D_sr_3) * (D_tv_sr / D_sr_3)));
                this.K_tv = (D_tv_sr / D_sr_3) * Math.sqrt(z_tv / (z_v * (m + 1)));
                System.out.println("\n-----------------------");
                System.out.println("z_tv = " + z_tv);
                System.out.println("K_tv = " + K_tv);
                System.out.println("-----------------------");
            }

            double u_tv_sr;
            // i — передаточное отношение редуктора, расположенного между роторами вентилятора и его турбины. Так как он отсутствует то i=1.
            double i_pered_chislo = 1;

            public void okryzhnaya_skorost_na_srednem_diametre_tyrbinbl() {
                this.u_tv_sr = Math.PI * D_tv_sr * (N_v / 60) * i_pered_chislo;
                System.out.println("\n-----------------------");
                System.out.println("i_pered_chislo = " + i_pered_chislo);
                System.out.println("u_tv_sr = " + u_tv_sr);
                System.out.println("-----------------------");
            }

            double Y_tv$;

            public void parametr_nagryzhennosti_stypeney_tyrbinbl_ventilyatora() {
                this.Y_tv$ = u_tv_sr * Math.sqrt((z_tv * Nu_tv$) / (2 * L_tv));
                System.out.println("\n-----------------------");
                System.out.println("Y_tv$ = " + Y_tv$);
                System.out.println("-----------------------");
            }

            double T_l;
            double k_sigma£;
            double k_sigma;
            double tau;
            double sigma_dl;
            double Rou_m;
            double sigma_p;
            double sigma_sum;
            double Fi = 0.68;

            private void prochnostnble_parametri_rabochey_lopatki_stypeni_tyrbinbl() {
                this.T_l = 0.95 * (T_tv$ + (u_tv_sr * u_tv_sr) / (2 * (k_gaza / (k_gaza - 1)) * R_gaza));
                //Зададим запас прочности рабочих лопаток k_sigma£ = sigma_dl / sigma_p >= 2
                this.k_sigma£ = 2;
                //материал ЭИ-4376
                this.tau = 400;
                this.sigma_dl = 290_000_000;
                this.Rou_m = 8470;
                this.sigma_p = sigma_dl / k_sigma£;
                this.sigma_sum = 8.2 * Rou_m * F_2_tv * Fi * N_v * i_pered_chislo;
                this.k_sigma = sigma_dl / sigma_sum;
                System.out.println("\n-----------------------");
                System.out.println("T_l = " + T_l);
                System.out.println("k_sigma£ = " + k_sigma£);
                System.out.println("k_sigma = " + k_sigma);
                System.out.println("tau = " + tau);
                System.out.println("sigma_dl = " + sigma_dl);
                System.out.println("Rou_m = " + Rou_m);
                System.out.println("sigma_p = " + sigma_p);
                System.out.println("sigma_sum = " + sigma_sum);
                System.out.println("-----------------------");
            }

            //Охлаждение не нужно тк температура корня лопатки ниже 1100К

            double D_sr_del_h;

            public void opredelenie_otn_vblsoti_lopatok_tyrbinbl() {
                this.D_sr_del_h = (2 * u_tv_sr * u_tv_sr * Rou_m * Fi) / sigma_p;
                System.out.println("\n-----------------------");
                System.out.println("D_sr_del_h = " + D_sr_del_h);
                System.out.println("-----------------------");
            }

            double D_2_tv_vt;

            public void vnytr_diametr_posled_styp_tyrbinbl() {
                this.D_2_tv_vt = D_tv_vblx_nar - 2 * h_2_tv;
                System.out.println("\n-----------------------");
                System.out.println("D_2_tv_vt = " + D_2_tv_vt);
                System.out.println("-----------------------");
            }

            public void calculate() {
                nar_diametr_tyrbinbl_na_vblxode();
                potrebnaya_vnytr_ydelnaya_rabota_tyrb_tyrboventilyatora();
                temperatyra_gaza_pered_tyrbinoy();
                calculate_poln_davleniy_v_tyrbine_ventilyatora();
                ploshad_kol_sechenia_na_vblxode_iz_tyrbinbl_ventilyatora();
                vblsota_lopatok_tyrbinbl_vent();
                sredniy_diameter_tyrbinbl_na_vblxode();
                ploshad_protochnoy_chasti_na_vxode_v_tyrbiny_ventilyatora();
                diametraln_razmeri_na_vxode_v_tyrbiny();
                chislo_stypeney_tyrbinbl();
                okryzhnaya_skorost_na_srednem_diametre_tyrbinbl();
                parametr_nagryzhennosti_stypeney_tyrbinbl_ventilyatora();
                prochnostnble_parametri_rabochey_lopatki_stypeni_tyrbinbl();
                opredelenie_otn_vblsoti_lopatok_tyrbinbl();
                vnytr_diametr_posled_styp_tyrbinbl();
            }

            class P4_Детальный_расчет_турбины_низкого_давления {

                double K_l = 0.05;
                double s_rk_1_oxl;

                public void opred_shirini_oxlazhdaemix_rab_reshetok() {
                    this.s_rk_1_oxl = K_l * D_sr_del_h * h_2_tv;
                    System.out.println("\n-----------------------");
                    System.out.println("K_l = " + K_l);
                    System.out.println("s_rk_1_oxl = " + s_rk_1_oxl);
                    System.out.println("-----------------------");
                }

                double s_ca_1_oxl;

                public void opred_shirini_soplov_apparatov() {
                    this.s_ca_1_oxl = 1.12 * s_rk_1_oxl;
                    System.out.println("\n-----------------------");
                    System.out.println("s_ca_1_oxl = " + s_ca_1_oxl);
                    System.out.println("-----------------------");
                }

                double s_rk_2_neoxl;
                double s_ca_2_neoxl;

                // Их ширина обычно составляет на 5…15% меньше, чем у охлаждаемых.
                public void opred_shirini_neoxlazhdaemix() {
                    this.s_rk_2_neoxl = s_rk_1_oxl / 1.07;
                    this.s_ca_2_neoxl = s_ca_1_oxl / 1.07;
                    System.out.println("\n-----------------------");
                    System.out.println("s_rk_2_neoxl = " + s_rk_2_neoxl);
                    System.out.println("s_ca_2_neoxl = " + s_ca_2_neoxl);
                    System.out.println("-----------------------");
                }

                double delta_s_1_oxl;
                double delta_s_2_neoxl;

                public void zazor_mezhdy_vench() {
                    this.delta_s_1_oxl = s_rk_1_oxl * 0.13;
                    this.delta_s_2_neoxl = s_rk_2_neoxl * 0.13;
                    System.out.println("\n-----------------------");
                    System.out.println("delta_s_1 = " + delta_s_1_oxl);
                    System.out.println("delta_s_2 = " + delta_s_2_neoxl);
                    System.out.println("-----------------------");
                }

                //Руководствуясь выше приведенными расчетами, можно заключить, что длина проточной части турбины составляет
                double l_t;

                //Можно расчет охлаждаемых и неохлаждаемых лопаток турбины
                public void dlina_protochn_chasti() {
                    this.l_t = s_rk_2_neoxl + s_ca_2_neoxl + 2 * delta_s_2_neoxl;
                    System.out.println("\n-----------------------");
                    System.out.println("l_t = " + l_t);
                    System.out.println("-----------------------");
                }

                //После этого, следует задаться рекомендованным значением угла выхода газа из сопловых аппаратов турбины.
                double alfa_1_st = 20;
                //Затем следует задаться скоростными коэффициентами на ступенях турбины.
                double fi = 0.98;
                double psi = 0.96;
                //Задаться степенью реактивности
                double Rou_st_tv = 0.35;

                double L_0_1_st;
                //— коэффициент нагруженности ступени.
                double Y_st$ = 0.5544;

                public void adiabat_rab_rashirenia() {
                    this.L_0_1_st = 1 / 2. * Math.pow(u_tv_sr / Y_st$, 2);
                    System.out.println("\n-----------------------");
                    System.out.println("alfa_1_st = " + alfa_1_st);
                    System.out.println("fi = " + fi);
                    System.out.println("psi = " + psi);
                    System.out.println("Rou_st_tv = " + Rou_st_tv);
                    System.out.println("L_0_1_st = " + L_0_1_st);
                    System.out.println("Y_st$ = " + Y_st$);
                    System.out.println("-----------------------");
                }

                double L_0_1;

                public void adiab_rabota_v_soplovom() {
                    this.L_0_1 = (1 - Rou_st_tv) * L_0_1_st;
                    System.out.println("\n-----------------------");
                    System.out.println("L_0_1 = " + L_0_1);
                    System.out.println("-----------------------");
                }

                double L_0_2;

                public void adiab_rabota_v_rabochem_kolese() {
                    this.L_0_2 = Rou_st_tv * L_0_1_st;
                    System.out.println("\n-----------------------");
                    System.out.println("L_0_2 = " + L_0_2);
                    System.out.println("-----------------------");
                }

                double c_1_t;

                public void teor_skorost_na_vblxode_iz_soplovogo() {
                    this.c_1_t = Math.sqrt(2 * L_0_1);
                    System.out.println("\n-----------------------");
                    System.out.println("c_1_t = " + c_1_t);
                    System.out.println("-----------------------");
                }

                double c_ad_t;

                public void ysl_skorosti_pri_adiab_rashrenii_gaza() {
                    this.c_ad_t = 1.415 * Math.sqrt(L_0_1_st);
                    System.out.println("\n-----------------------");
                    System.out.println("c_ad_t = " + c_ad_t);
                    System.out.println("-----------------------");
                }

                double lambda_ad;

                public void prived_skorost() {
                    this.lambda_ad = c_ad_t / (Math.sqrt(2 * k_gaza / (k_gaza + 1) * R_gaza * T_vblx_tk$));
                    System.out.println("\n-----------------------");
                    System.out.println("lambda_ad = " + lambda_ad);
                    System.out.println("-----------------------");
                }

                double pi_lambda_ad = 0.6681;
                double P_2_t;

                public void staticheskogo_davlenia_za_tyrbinoi() {
                    this.P_2_t = P_tk$ * pi_lambda_ad;
                    System.out.println("\n-----------------------");
                    System.out.println("P_2_t = " + P_2_t);
                    System.out.println("pi_lambda_ad = " + pi_lambda_ad);
                    System.out.println("-----------------------");
                }

                double c_1_T;

                public void deistvit_skorost_na_vblxode_iz_sopla() {
                    this.c_1_T = fi * c_1_t;
                    System.out.println("\n-----------------------");
                    System.out.println("c_1_T = " + c_1_T);
                    System.out.println("-----------------------");
                }

                double T_1;

                private void stat_temp_gaza_za_sopl_tyrb() {
                    this.T_1 = T_vblx_tk$ - (fi * fi * L_0_1) / (k_gaza / (k_gaza - 1) * R_gaza);
                    System.out.println("\n-----------------------");
                    System.out.println("T_1 = " + T_1);
                    System.out.println("-----------------------");
                }

                double lambda_1_t;

                private void teor_prived_skorosti_na_vblxode_iz_sopl_apparata() {
                    this.lambda_1_t = c_1_t / Math.sqrt(2 * k_gaza / (k_gaza + 1) * R_gaza * T_vblx_tk$);
                    System.out.println("\n-----------------------");
                    System.out.println("lambda_1_t = " + lambda_1_t);
                    System.out.println("-----------------------");
                }

                //Зная приведенную скорость на выходе из соплового аппарата, можно определить соответствующие им газодинамические функции
                //lambda_1_t = 0.6652956049660866
                double pi_lambda_1_t = 0.7737;
                double q_lambda_1_t = 0.8645;

                double P_1_tv;

                public void stat_davlen_za_soplov_app_stypeni_tyrbinbl() {
                    this.P_1_tv = pi_lambda_1_t * P_tk$;
                    System.out.println("\n-----------------------");
                    System.out.println("P_1_tv = " + P_1_tv);
                    System.out.println("pi_lambda_1_t = " + pi_lambda_1_t);
                    System.out.println("q_lambda_1_t = " + q_lambda_1_t);
                    System.out.println("lambda_1_t = " + lambda_1_t);
                    System.out.println("-----------------------");
                }

                double Rou_1;

                public void opred_plotn_gaza() {
                    this.Rou_1 = P_1_tv / (R_gaza * T_1);
                    System.out.println("\n-----------------------");
                    System.out.println("Rou_1 = " + Rou_1);
                    System.out.println("-----------------------");
                }

                double sin_alfa_1_tv;
                double alfa_1_tv;
                double h_1 = 0.0548;

                private void ygl_vblxoda_iz_soplovogo_apparata() {
                    this.sin_alfa_1_tv = G_t / (Math.PI * D_tv_sr * h_1 * c_1_t * Rou_1);
                    this.alfa_1_tv = fromRadToGrad(Math.asin(sin_alfa_1_tv));
                    System.out.println("\n-----------------------");
                    System.out.println("h_1 = " + h_1);
                    System.out.println("sin_alfa_1_tv = " + sin_alfa_1_tv);
                    System.out.println("alfa_1_tv = " + alfa_1_tv);
                    System.out.println("-----------------------");
                }

                double Rou_st_k;

                public void stepen_reaktivnosti_y_kornya_lopatki() {
                    this.Rou_st_k = 1 - (1 - Rou_st_tv)
                            * (Math.pow(D_tv_sr / (D_tv_sr - h_1), 2)
                            * Math.pow(Math.cos(fromGradToRad(alfa_1_tv)), 2)
                            + Math.pow(Math.sin(fromGradToRad(alfa_1_tv)), 2));
                    System.out.println("\n-----------------------");
                    System.out.println("Rou_st_k = " + Rou_st_k);
                    System.out.println("cos(alfa_1_tv)^2 = " + Math.pow(Math.cos(fromGradToRad(alfa_1_tv)), 2));
                    System.out.println("sin(alfa_1_tv)^2 = " + Math.pow(Math.sin(fromGradToRad(alfa_1_tv)), 2));
                    System.out.println("-----------------------");
                }

                double W_1_T;
                double cos_alfa_1_tv;

                public void opredelenie_skorosti_potoka__gaza_na_vxode_v_rabochee_koleso() {
                    this.cos_alfa_1_tv = Math.cos(fromGradToRad(alfa_1_tv));
                    this.W_1_T = Math.sqrt(c_1_T * c_1_T + u_tv_sr * u_tv_sr - 2 * c_1_T * u_tv_sr * cos_alfa_1_tv);
                    System.out.println("\n-----------------------");
                    System.out.println("W_1_T = " + W_1_T);
                    System.out.println("cos_alfa_1_tv = " + cos_alfa_1_tv);
                    System.out.println("-----------------------");
                }

                double sin_betta_1_tv;
                double betta_1_tv;

                public void ygla_vxoda_potoka_na_rabochyi_reshetky_tyrbinbl() {
                    this.sin_betta_1_tv = (c_1_T * sin_alfa_1_tv) / W_1_T;
                    this.betta_1_tv = fromRadToGrad(Math.asin(sin_betta_1_tv));
                    System.out.println("\n-----------------------");
                    System.out.println("sin_betta_1_tv = " + sin_betta_1_tv);
                    System.out.println("betta_1_tv = " + betta_1_tv);
                    System.out.println("-----------------------");
                }

                double W_2_T;

                public void skorosti_gaza_na_vblxode_iz_rabochey_reshetki_turbinbl() {
                    this.W_2_T = psi * Math.sqrt(W_1_T * W_1_T + 2 * L_0_2);
                    System.out.println("\n-----------------------");
                    System.out.println("W_2_T = " + W_2_T);
                    System.out.println("-----------------------");
                }

                double T_w_1_T$;

                public void temp_tormozheniya_na_vblxode_v_otn_dvizhenii_na_vxode_v_reshetky() {
                    this.T_w_1_T$ = T_vblx_tk$ - (c_1_T * c_1_T - W_1_T * W_1_T) / (2 * k_gaza / (k_gaza - 1) * R_gaza);
                    System.out.println("\n-----------------------");
                    System.out.println("T_w_1_T$ = " + T_w_1_T$);
                    System.out.println("-----------------------");
                }

                double lambda_w_2_T;

                public void otnos_skorost_na_vblxode_iz_tyrbinbl() {
                    this.lambda_w_2_T = W_2_T / Math.sqrt(2 * k_gaza / (k_gaza + 1) * R_gaza * T_w_1_T$);
                    System.out.println("\n-----------------------");
                    System.out.println("lambda_w_2_T = " + lambda_w_2_T);
                    System.out.println("-----------------------");
                }

                //lambda_w_2_T = 0.225  -
                //double q_lambda_w_2_T = 0.3391;
                //double pi_lambda_w_2_T = 0.8478;

                //lambda_w_2_T = 0.599  +
                double q_lambda_w_2_T = 0.8133;
                double pi_lambda_w_2_T = 0.8098;

                double P_w_2_T$;

                public void polnoe_davlenie_v_otnosit_dvizhenii() {
                    this.P_w_2_T$ = P_2_t / pi_lambda_w_2_T;
                    System.out.println("\n-----------------------");
                    System.out.println("P_w_2_T$ = " + P_w_2_T$);
                    System.out.println("-----------------------");
                }

                double sin_betta_2_tv;
                //по чертежу - высота лопаток на выходе из рабочего колеса турбины
                double h_2 = 0.069;
                double betta_2_tv;

                public void opredelenie_ygla_vblxoda_iz_rab_reshetki_tyrb_v_otn_dvizh() {
                    this.sin_betta_2_tv = G_t * Math.sqrt(T_w_1_T$) / (Math.PI * D_tv_sr * h_2 * P_w_2_T$ * q_lambda_w_2_T * s_g);
                    this.betta_2_tv = fromRadToGrad(Math.asin(sin_betta_2_tv));
                    System.out.println("\n-----------------------");
                    System.out.println("h_2 = " + h_2);
                    System.out.println("sin_betta_2_tv = " + sin_betta_2_tv);
                    System.out.println("betta_2_tv = " + betta_2_tv);
                    System.out.println("-----------------------");
                }

                double c_2_T;
                double cos_betta_2_tv;

                public void skorost_potoka_za_rab_kolesom_tyrbinbl() {
                    this.cos_betta_2_tv = Math.cos(fromGradToRad(betta_2_tv));
                    this.c_2_T = Math.sqrt(W_2_T * W_2_T + u_tv_sr * u_tv_sr - 2 * W_2_T * u_tv_sr * cos_betta_2_tv);
                    System.out.println("\n-----------------------");
                    System.out.println("cos_betta_2_tv = " + cos_betta_2_tv);
                    System.out.println("c_2_T = " + c_2_T);
                    System.out.println("-----------------------");
                }

                double sin_alfa_2_tv;
                double alfa_2_tv;

                public void ygl_absol_skorosti_potoka_za_rab_kolesom() {
                    this.sin_alfa_2_tv = W_2_T * sin_betta_2_tv / c_2_T;
                    this.alfa_2_tv = fromRadToGrad(Math.asin(sin_alfa_2_tv));
                    System.out.println("\n-----------------------");
                    System.out.println("sin_alfa_2_tv = " + sin_alfa_2_tv);
                    System.out.println("alfa_2_tv = " + alfa_2_tv);
                    System.out.println("-----------------------");
                }

                double delta_t_opt_ca;
                double K_ca;
                double sin_alfa_0_tv = 1;
                double alfa_0_tv = 90;
                double delta_t_opt_rk;
                double K_rk;


                public void otnos_optimal_shag_reshetki() {
                    this.delta_t_opt_ca = -0.625 * lambda_1_t * lambda_1_t + 0.48 * lambda_1_t + 0.016;
                    this.K_ca = sin_alfa_0_tv / sin_alfa_1_tv;
                    this.delta_t_opt_rk = -0.625 * lambda_w_2_T * lambda_w_2_T + 0.48 * lambda_w_2_T + 0.016;
                    this.K_rk = sin_betta_1_tv / sin_betta_2_tv;
                    System.out.println("\n-----------------------");
                    System.out.println("delta_t_opt_ca = " + delta_t_opt_ca);
                    System.out.println("K_ca = " + K_ca);
                    System.out.println("delta_t_opt_rk = " + delta_t_opt_rk);
                    System.out.println("K_rk = " + K_rk);
                    System.out.println("-----------------------");
                }

                //ОТ К ЗАВИСЯТ ДАЛЬНЕЙШИЕ ФОРМУЛЫ
                // ТИП ЛОПАТОК ТОЖЕ ВЛИЯЕТ НА ФОРМУЛЫ ОХЛ ИЛИ НЕТ

                double tetta_rk;
                double tetta_ca;

                public void ygol_povorota_v_reshetkax() {
                    this.tetta_rk = (180 - (betta_1_tv + betta_2_tv)) * Math.PI / 180;
                    this.tetta_ca = (180 - (alfa_1_tv + alfa_2_tv)) * Math.PI / 180;
                    System.out.println("\n-----------------------");
                    System.out.println("tetta_rk = " + tetta_rk);
                    System.out.println("tetta_ca = " + tetta_ca);
                    System.out.println("-----------------------");
                }

                double t_opt_0_ca;
                double t_opt_0_rk;
                double Y_yst_ca;
                double Y_yst_rk;
                double b_ca;
                double b_rk;
                // 0.007 - 0.015
                double koef_r_vbx_ca = 0.015;
                double r_vbx_ca;
                double s_vbx_ca;
                // 0.015 - 0.02
                double koef_r_vbx_rk = 0.015;
                double r_vbx_rk;
                double s_vbx_rk;
                double K_kr_ca;
                double K_kr_rk;
                double t_opt_ca;
                double t_opt_rk;
                double sin_Y_yst_ca;
                double sin_Y_yst_rk;

                public void optim_shag() {
                    this.t_opt_0_ca = (0.327 / (Math.pow(K_ca, 0.371) * Math.pow(tetta_ca, 1 / 3.0))) - 0.994 / Math.pow(K_ca, 0.385) + 1.314;
                    this.t_opt_0_rk = (0.327 / (Math.pow(K_rk, 0.371) * Math.pow(tetta_rk, 1 / 3.0))) - 0.994 / Math.pow(K_rk, 0.385) + 1.314;
                    this.Y_yst_ca = 68.7 + 0.000933 * (alfa_1_tv - alfa_2_tv) - 0.006052 * (alfa_1_tv - alfa_2_tv) * (alfa_1_tv - alfa_2_tv);
                    this.Y_yst_rk = 68.7 + 0.000933 * (betta_1_tv - betta_2_tv) - 0.006052 * (betta_1_tv - betta_2_tv) * (betta_1_tv - betta_2_tv);
                    this.sin_Y_yst_ca = Math.sin(fromGradToRad(Y_yst_ca));
                    this.sin_Y_yst_rk = Math.sin(fromGradToRad(Y_yst_rk));
                    this.b_ca = s_ca_2_neoxl / sin_Y_yst_ca;
                    this.b_rk = s_rk_2_neoxl / sin_Y_yst_rk;
                    this.r_vbx_ca = koef_r_vbx_ca * b_ca;
                    this.r_vbx_rk = koef_r_vbx_rk * b_rk;
                    this.s_vbx_ca = 2 * r_vbx_ca / b_ca;
                    this.s_vbx_rk = 2 * r_vbx_rk / b_rk;
                    this.K_kr_ca = 1 - 1.5 * s_vbx_ca * s_vbx_ca + (3.75 - t_opt_0_ca - 0.6) * s_vbx_ca;
                    this.K_kr_rk = 1 - 1.5 * s_vbx_rk * s_vbx_rk + (3.75 - t_opt_0_rk - 0.6) * s_vbx_rk;
                    this.t_opt_ca = (1 - delta_t_opt_ca) * K_kr_ca * t_opt_0_ca;
                    this.t_opt_rk = (1 + delta_t_opt_rk) * K_kr_rk * t_opt_0_rk;
                    System.out.println("\n-----------------------");
                    System.out.println("t_opt_0_ca = " + t_opt_0_ca);
                    System.out.println("t_opt_0_rk = " + t_opt_0_rk);
                    System.out.println("Y_yst_ca = " + Y_yst_ca);
                    System.out.println("Y_yst_rk = " + Y_yst_rk);
                    System.out.println("b_ca = " + b_ca);
                    System.out.println("b_rk = " + b_rk);
                    System.out.println("koef_r_vbx_ca = " + koef_r_vbx_ca);
                    System.out.println("r_vbx_ca = " + r_vbx_ca);
                    System.out.println("koef_r_vbx_rk = " + koef_r_vbx_rk);
                    System.out.println("r_vbx_rk = " + r_vbx_rk);
                    System.out.println("s_vbx_ca = " + s_vbx_ca);
                    System.out.println("s_vbx_rk = " + s_vbx_rk);
                    System.out.println("K_kr_ca = " + K_kr_ca);
                    System.out.println("K_kr_rk = " + K_kr_rk);
                    System.out.println("t_opt_ca = " + t_opt_ca);
                    System.out.println("t_opt_rk = " + t_opt_rk);
                    System.out.println("sin_Y_yst_ca = " + sin_Y_yst_ca);
                    System.out.println("sin_Y_yst_rk = " + sin_Y_yst_rk);
                    System.out.println("-----------------------");
                }

                double t_ca;
                double t_rk;
                int z_ca;
                int z_rk;

                public void shag_reshetki() {
                    System.out.println("\n-----------------------");
                    this.t_ca = t_opt_ca * b_ca;
                    System.out.println("t_ca [ДО]= " + t_ca);
                    this.t_rk = t_opt_rk * b_rk;
                    System.out.println("t_rk [ДО] = " + t_rk);
                    this.z_ca = (int) (Math.PI * D_tv_sr / t_ca);
                    this.t_ca = (Math.PI * D_tv_sr / z_ca);
                    System.out.println("t_ca [ПОСЛЕ] = " + t_ca);
                    this.z_rk = (int) (Math.PI * D_tv_sr / t_rk);
                    this.t_rk = (Math.PI * D_tv_sr / z_rk);
                    System.out.println("t_rk [ПОСЛЕ] = " + t_rk);
                    System.out.println("z_ca = " + z_ca);
                    System.out.println("z_rk = " + z_rk);
                    System.out.println("-----------------------");
                }

                public void calculate() {
                    opred_shirini_oxlazhdaemix_rab_reshetok();
                    opred_shirini_soplov_apparatov();
                    opred_shirini_neoxlazhdaemix();

                    ////////////////////
                    zazor_mezhdy_vench();
                    dlina_protochn_chasti();
                    adiabat_rab_rashirenia();
                    adiab_rabota_v_soplovom();
                    adiab_rabota_v_rabochem_kolese();
                    teor_skorost_na_vblxode_iz_soplovogo();
                    ysl_skorosti_pri_adiab_rashrenii_gaza();
                    prived_skorost();
                    staticheskogo_davlenia_za_tyrbinoi();
                    deistvit_skorost_na_vblxode_iz_sopla();
                    stat_temp_gaza_za_sopl_tyrb();
                    teor_prived_skorosti_na_vblxode_iz_sopl_apparata();
                    stat_davlen_za_soplov_app_stypeni_tyrbinbl();
                    opred_plotn_gaza();
                    ygl_vblxoda_iz_soplovogo_apparata();
                    stepen_reaktivnosti_y_kornya_lopatki();
                    opredelenie_skorosti_potoka__gaza_na_vxode_v_rabochee_koleso();
                    ygla_vxoda_potoka_na_rabochyi_reshetky_tyrbinbl();
                    skorosti_gaza_na_vblxode_iz_rabochey_reshetki_turbinbl();
                    temp_tormozheniya_na_vblxode_v_otn_dvizhenii_na_vxode_v_reshetky();
                    otnos_skorost_na_vblxode_iz_tyrbinbl();
                    polnoe_davlenie_v_otnosit_dvizhenii();
                    opredelenie_ygla_vblxoda_iz_rab_reshetki_tyrb_v_otn_dvizh();
                    skorost_potoka_za_rab_kolesom_tyrbinbl();
                    ygl_absol_skorosti_potoka_za_rab_kolesom();
                    otnos_optimal_shag_reshetki();
                    ygol_povorota_v_reshetkax();
                    optim_shag();
                    shag_reshetki();
                }

                class Р5_Расчет_параметров_потока_по_радиусу_ступени_турбины {

                    //cюда формулу для нахождения, слайдеры
                    // +- (2-4)мм корень
                    double r_1 = 0.294;
                    //Dcp / 2 средний
                    double r_2 = 0.31075;
                    // +- (2-4)мм профиль
                    double r_3 = 0.328;

                    //выбор закона закрутки

                    //среднее взято как входное и выходное на первой ступени турбины
                    // c_1_T
                    double c_1_a_sr = 423.194;

                    // c_2_T
                    double c_2_a_sr = 212.143;

                    //Так как закон закрутки с постоянной циркуляцией вектора скорости, то
                    double c_1_a_1;
                    double c_1_a_2;
                    double c_1_a_3;
                    double c_2_a_1;
                    double c_2_a_2;
                    double c_2_a_3;

                    public void zakon_zakrytki() {
                        this.c_1_a_1 = c_1_a_sr;
                        this.c_1_a_2 = c_1_a_sr;
                        this.c_1_a_3 = c_1_a_sr;
                        this.c_2_a_1 = c_2_a_sr;
                        this.c_2_a_2 = c_2_a_sr;
                        this.c_2_a_3 = c_2_a_sr;
                        System.out.println("\n-----------------------");
                        System.out.println("c_1_a_1 = " + c_1_a_1);
                        System.out.println("c_1_a_2 = " + c_1_a_2);
                        System.out.println("c_1_a_3 = " + c_1_a_3);
                        System.out.println("c_2_a_1 = " + c_2_a_1);
                        System.out.println("c_2_a_2 = " + c_2_a_2);
                        System.out.println("c_2_a_3 = " + c_2_a_3);
                        System.out.println("-----------------------");
                    }

                    //переферийный радиус
                    double r_k = 0.31075 * 2;

                    double r_otn_sr;
                    double r_otn_1;
                    double r_otn_2;
                    double r_otn_3;

                    public void otnosit_radiysi() {
                        this.r_otn_sr = r_2 / r_k;
                        this.r_otn_1 = r_1 / r_k;
                        this.r_otn_2 = r_2 / r_k;
                        this.r_otn_3 = r_3 / r_k;

                        System.out.println("\n-----------------------");
                        System.out.println("r_k = " + r_k);
                        System.out.println("r_otn_sr = " + r_otn_sr);
                        System.out.println("r_otn_1 = " + r_otn_1);
                        System.out.println("r_otn_2 = " + r_otn_2);
                        System.out.println("r_otn_3 = " + r_otn_3);
                        System.out.println("-----------------------");
                    }

                    // они задаются что и определяет закон закрутки
                    double c_1_u_sr;
                    double c_2_u_sr;

                    double c_1_u_1;
                    double c_1_u_2;
                    double c_1_u_3;
                    double c_2_u_1;
                    double c_2_u_2;
                    double c_2_u_3;

                    public void opredel_okryzhnoi_sostavl_absol_skor_vdol_radiysa_na_vxode() {
                        this.c_1_u_sr = u_tv_sr;
                        this.c_2_u_sr = u_tv_sr;

                      /*  this.c_1_u_1 = c_1_u_s
                        this.c_2_u_1
                        this.c_1_u_2
                        this.c_2_u_2
                        this.c_1_u_3
                        this.c_2_u_3*/
                    }

                    public void calculate() {
                    }
                }
            }
        }
    }
}
