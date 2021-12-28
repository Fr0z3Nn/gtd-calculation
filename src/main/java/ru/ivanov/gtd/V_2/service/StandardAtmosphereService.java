package ru.ivanov.gtd.V_2.service;

import org.springframework.stereotype.Service;
import ru.ivanov.gtd.V_2.domain.StandardAtmosphereDto;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author Sergey Ivanov
 * created on 08.11.2021
 */
@Service
public class StandardAtmosphereService {
    public StandardAtmosphereDto processStandardAtmosphere(StandardAtmosphereDto standardAtmosphereDto) {

        boolean isNotExistH_n = standardAtmosphereDto.getH_n() == null;

        boolean isNotExistT_vx$Or_P_n$ = standardAtmosphereDto.getT_vx$() == null
                || standardAtmosphereDto.getP_n$() == null;


        if (isNotExistH_n && isNotExistT_vx$Or_P_n$) {
            throw new RuntimeException("Должны быть заданы начальные параметры стандартной атмосферы");
        }
        if (isNotExistH_n) {
            return standardAtmosphereDto;
        } else {
            return calculateStandardParamByH_n(standardAtmosphereDto.getH_n());
        }
    }

    private StandardAtmosphereDto calculateStandardParamByH_n(double H_n) {
        Integer key = (int) H_n;
        List<Double> standardParam = Arrays.stream(standardAtmosphereMap.get(key).split("-"))
                .map(Double::valueOf)
                .collect(Collectors.toList());
        return StandardAtmosphereDto.builder()
                .H_n(H_n)
                .T_vx$(standardParam.get(0))
                .P_n$(standardParam.get(1))
                .build();
    }

    //ЗДЕСЬ ПОЛНЫЙ СПИСОК ИЛИ ПОЛИНОМ
    private Map<Integer, String> standardAtmosphereMap =
            Map.of(0, "288-101325");
}
