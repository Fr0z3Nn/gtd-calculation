package ru.ivanov.gtd.V_2.service;

import lombok.RequiredArgsConstructor;
import org.springframework.stereotype.Service;
import ru.ivanov.gtd.V_2.domain.fan.FanRequestDto;

/**
 * @author Sergey Ivanov
 * created on 08.11.2021
 */
@Service
@RequiredArgsConstructor
public class PreliminaryFanCalculationService {

    private final StandardAtmosphereService atmosphereStandardService;

    public String calculate(FanRequestDto fanRequestDto){
       // ObjectMapper objectMapper = new ObjectMapper();

        return atmosphereStandardService.processStandardAtmosphere(fanRequestDto.getStandardAtmosphere()).toString();
    }

}
