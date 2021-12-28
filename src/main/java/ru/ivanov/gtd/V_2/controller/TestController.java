package ru.ivanov.gtd.V_2.controller;

import lombok.RequiredArgsConstructor;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RestController;
import ru.ivanov.gtd.V_2.domain.fan.FanRequestDto;
import ru.ivanov.gtd.V_2.service.PreliminaryFanCalculationService;

import javax.validation.Valid;

/**
 * @author Sergey Ivanov
 * created on 08.11.2021
 */
@RestController
@RequiredArgsConstructor
public class TestController {

    private final PreliminaryFanCalculationService preliminaryFanCalculationService;

    @GetMapping("/test")
    public String test(@RequestBody @Valid FanRequestDto fanRequestDto) {
        System.out.println(fanRequestDto.toString());
        return preliminaryFanCalculationService.calculate(fanRequestDto);
    }
}
