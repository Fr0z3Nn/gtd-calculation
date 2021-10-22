package ru.ivanov.gtd.Ship;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.Configuration;
import org.springframework.web.servlet.config.annotation.EnableWebMvc;
import org.springframework.web.servlet.config.annotation.WebMvcConfigurer;

/**
 * @author Sergey Ivanov
 * created on 20.10.2021
 */
@Configuration
@EnableWebMvc
@ComponentScan(basePackages = "ru.ivanov.gtd")
public class AppConfig implements WebMvcConfigurer {
    @Bean
    public InjectValueAnnotationPostBeanProcessor injectValueAnnotationPostBeanProcessor(){
        return new InjectValueAnnotationPostBeanProcessor();
    }

    @Bean
    public ReplaceAnnotationPostBeanProcessor replaceAnnotationPostBeanProcessor(){
        return new ReplaceAnnotationPostBeanProcessor();
    }
}
