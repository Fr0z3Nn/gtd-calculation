package ru.ivanov.gtd.V_2.domain.fan;

import com.fasterxml.jackson.annotation.JsonProperty;
import lombok.*;
import ru.ivanov.gtd.V_2.domain.StandardAtmosphereDto;

import javax.validation.constraints.DecimalMax;
import javax.validation.constraints.DecimalMin;
import javax.validation.constraints.NotNull;

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
public class FanRequestDto {

    @DecimalMin("0.6")
    @DecimalMax("0.75")
    private double lambda_vx;

    @DecimalMin("0.3")
    @DecimalMax("0.45")
    private double otn_d_vt;

    private double G_v_sum;

    private double m;

    private double Pi_sum$;

    private double Pi_v$;

    @NotNull
    private StandardAtmosphereDto standardAtmosphere;
}
