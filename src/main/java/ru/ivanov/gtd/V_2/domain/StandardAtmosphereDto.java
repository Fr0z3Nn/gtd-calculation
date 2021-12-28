package ru.ivanov.gtd.V_2.domain;

import lombok.*;

/**
 * @author Sergey Ivanov
 * created on 08.11.2021
 */
@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@Builder
@ToString
public class StandardAtmosphereDto {

    private Double H_n;

    private Double T_vx$;

    private Double P_n$;
}
