package ru.ivanov.gtd;

import org.springframework.boot.autoconfigure.SpringBootApplication;
import ru.ivanov.gtd.Турбовентилятор.P1_Предварительный_расчет_вентилятора;

import java.util.Optional;
import java.util.Random;


@SpringBootApplication

public class GtdCalculationApplication {

    public static void main(String[] args) {
        var p1 = new P1_Предварительный_расчет_вентилятора();
        p1.calculate();
        //SpringApplication.run(GtdCalculationApplication.class, args);
    }

}




