package ru.ivanov.gtd.Ship;

import org.springframework.stereotype.Component;

import javax.annotation.PostConstruct;

/**
 * @author Sergey Ivanov
 * created on 20.10.2021
 */
@Component
public class Manager {

    @InjectValue("Viktor")
    @Replace("workPlace")
    private String name;

    @InjectValue("PAO")
    @Replace("name")
    private String workPlace;

    public Manager() {
        System.out.println("Объект создан");
    }

    @PostConstruct
    private void showCurrentValues(){
        System.out.println("Объект инициализирован" + this.toString());
    }

    @Override
    public String toString() {
        return "Manager{" +
                "name='" + name + '\'' +
                ", workPlace='" + workPlace + '\'' +
                '}';
    }
}
