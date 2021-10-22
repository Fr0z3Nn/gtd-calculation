package ru.ivanov.gtd.Ship;

import org.springframework.beans.factory.config.BeanPostProcessor;
import org.springframework.util.ReflectionUtils;

import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Sergey Ivanov
 * created on 20.10.2021
 */

public class ReplaceAnnotationPostBeanProcessor implements BeanPostProcessor {

    Map<Field, String> originalNames = new HashMap<>();

    public ReplaceAnnotationPostBeanProcessor() {
    }

    @Override
    public Object postProcessAfterInitialization(Object bean, String beanName) {
        Class<?> aClass = bean.getClass();
        String fieldValue;
        Field[] declaredFields = aClass.getDeclaredFields();
        for (Field field : declaredFields) {
            Replace annotation = field.getAnnotation(Replace.class);
            if (annotation != null) {
                String annotationValue = annotation.value();
                try {
                    field.setAccessible(true);
                    Field declaredField = aClass.getDeclaredField(annotationValue);
                    declaredField.setAccessible(true);
                    fieldValue = (String) ReflectionUtils.getField(declaredField, bean);
                    String previousValue = (String) ReflectionUtils.getField(field, bean);
                    String originalNamesOrDefault = originalNames.getOrDefault(declaredField, fieldValue);
                    ReflectionUtils.setField(field, bean, originalNamesOrDefault);
                    System.out.println("Изменено значение поля " + field.getName() + " на " + originalNamesOrDefault);
                    originalNames.put(field, previousValue);

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        return bean;
    }
}
