package ru.ivanov.gtd.Турбовентилятор;

import java.util.stream.IntStream;

import static ru.ivanov.gtd.Util.*;
import static ru.ivanov.gtd.Турбовентилятор.P00_Исходные_данные.*;
import static ru.ivanov.gtd.Турбовентилятор.P0_Константы.*;

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
        this.L_v = (k / (k - 1)) * R * T_n$ * (Math.pow(P00_Исходные_данные.Pi_v$, (k - 1) / k) - 1) * (1 / Nu_v$);
        this.T_k_vx$ = T_vx$ + (L_v / ((k / (k - 1)) * R));
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
        this.T_1_2$ = T_1_1$ + (L_v_1 / (k / (k - 1) * R));
        this.T_2_1$ = T_1_2$;
        this.T_3_1$ = T_1_2$;

        this.T_1_3$ = T_1_2$ + (L_v_2 / (k / (k - 1) * R));
        this.T_2_2$ = T_1_3$;
        this.T_3_2$ = T_1_3$;

        this.T_2_3$ = T_1_3$ + (L_v_3 / (k / (k - 1) * R));
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
        this.Pi_v_1$ = Math.pow((L_v_1 * Nu_v_1$) / ((k / (k - 1)) * R * T_1_1$) + 1, (k / (k - 1)));
        this.Pi_v_2$ = Math.pow((L_v_2 * Nu_v_1$) / ((k / (k - 1)) * R * T_1_2$) + 1, (k / (k - 1)));
        this.Pi_v_3$ = Math.pow((L_v_3 * Nu_v_1$) / ((k / (k - 1)) * R * T_1_3$) + 1, (k / (k - 1)));
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
        this.alfa_1_kr_1 = Math.sqrt(2 * k / (k + 1) * R * T_1_1$);
        this.alfa_1_kr_2 = Math.sqrt(2 * k / (k + 1) * R * T_1_2$);
        this.alfa_1_kr_3 = Math.sqrt(2 * k / (k + 1) * R * T_1_3$);
        this.alfa_1_kr_vblx = Math.sqrt(2 * k / (k + 1) * R * T_vblx$);
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
        this.alfa_1_kr_kompressor = Math.sqrt(2 * k / (k + 1) * R * T_k_vx$);
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
        double alfa_1_1 = Math.sqrt(((k + 1) / 2) * tau_lamda_1_1);

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
        this.lamda_1_1_f = c_1_1 / Math.sqrt(2 * k / (k + 1) * R * T_1_1$);
        this.lamda_1_2_f = c_1_2 / Math.sqrt(2 * k / (k + 1) * R * T_1_2$);
        this.lamda_1_3_f = c_1_3 / Math.sqrt(2 * k / (k + 1) * R * T_1_3$);
        this.lamda_1_vblx_f = c_1_3 / Math.sqrt(2 * k / (k + 1) * R * T_vblx$);
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
        var predVaritTND = new Р3_Предварительный_расчет_турбины_низкого_давления();
        predVaritTND.calculate();
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
            this.a_1_1 = alfa_1_kr_1 * Math.sqrt((k + 1) / 2 * tau_lambda_1_1);
            this.a_1_2 = alfa_1_kr_2 * Math.sqrt((k + 1) / 2 * tau_lambda_1_2);
            this.a_1_3 = alfa_1_kr_3 * Math.sqrt((k + 1) / 2 * tau_lambda_1_3);
            this.a_1_vblx = alfa_1_kr_vblx * Math.sqrt((k + 1) / 2 * tau_lambda_v_vblx);
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
            this.a_2_1 = alfa_1_kr_2 * Math.sqrt((k + 1) / 2 * tau_lambda_2_1);
            this.a_2_2 = alfa_1_kr_3 * Math.sqrt((k + 1) / 2 * tau_lambda_2_2);
            this.a_2_3 = alfa_1_kr_vblx * Math.sqrt((k + 1) / 2 * tau_lambda_2_3);
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
            this.P_2_1$ = P_1_1$ * Math.pow((1 + ((L_v_1 * Nu_rk$) / (k / (k - 1) * R * T_vx$))), k / (k - 1));
            this.P_2_2$ = P_1_2$ * Math.pow((1 + ((L_v_2 * Nu_rk$) / (k / (k - 1) * R * T_1_2$))), k / (k - 1));
            this.P_2_3$ = P_1_3$ * Math.pow((1 + ((L_v_3 * Nu_rk$) / (k / (k - 1) * R * T_1_3$))), k / (k - 1));
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
        public void opredelenie_chisla_lopatok_spryamlyash_apparatov(){
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

        public void opred_hord_streml_apparatov(){
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
    }

    class Р3_Предварительный_расчет_турбины_низкого_давления{

        public void calculate(){}
    }
}
