package ru.ivanov.gtd.V_1;

/**
 * @author Sergey Ivanov
 * created on 07.10.2021
 */
public class Util {
    public static double fromGradToRad(double grad) {
        return (grad * Math.PI) / 180;
    }

    public static double fromRadToGrad(double grad) {
        return (grad * 180) / Math.PI;
    }

    //gradysi
    public static double arcctg(double x){
        return  fromRadToGrad((Math.PI / 2) - Math.atan(x));
       //return  fromRadToGrad((Math.PI / 2) - Math.atan(x));
    }

    //x в градусах ответ
    public static double ctg(double x){
        return 1 / Math.tan(fromGradToRad(x));
    }
}
