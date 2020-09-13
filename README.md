# Осциллятор Дуффинга. Реализация на Python

Осциллятор Дуффинга  — простейшая одномерная нелинейная система. Представляет собой одномерную частицу, движущуюся в потенциале. 
При b = 0 система сводится к обычному линейному осциллятору. Особенностью осциллятора Дуффинга является возможность получения хаотической динамики.

Подробнее: https://ru.wikipedia.org/wiki/%D0%9E%D1%81%D1%86%D0%B8%D0%BB%D0%BB%D1%8F%D1%82%D0%BE%D1%80_%D0%94%D1%83%D1%84%D1%84%D0%B8%D0%BD%D0%B3%D0%B0

Этот код имитирует осциллятор затухания:
Осциллятор с затухающими гармониками в двухъямном потенциале.

  F = -Gamma * ( dx / dt ) + 2 * a * x - 4 * b * ( x ^ 3 ) + F_0 * cos( Omega * t )

Нелинейное дифференциальное уравнение второго порядка, решенное численно с помощью разложения Тейлора.

Для текущего набора параметров движение хаотическое, т.е. движение сильно зависит от начальных условий. 
Дополнительно фиксированного периода движения не наблюдается. Сюжет Пуанкаре - фрактал.

# Пример работы программы с стадартными параметрами
![alt tag](https://github.com/PC-SET/Duffing-Oscillator/blob/master/Example.png?raw=true "Пример")​

# Пример работы программы с хаотическим поведением
![alt tag](https://github.com/PC-SET/Duffing-Oscillator/blob/master/Example2.png?raw=true "Пример")​

