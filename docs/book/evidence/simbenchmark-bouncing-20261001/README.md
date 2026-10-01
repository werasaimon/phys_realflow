# SimBenchmark Bouncing: срез 2026-10-01

Параметры и математика — в [главе 13](../../13-external-benchmark.md).
Это адаптация одного авторского опыта, не прохождение всего SimBenchmark.
Все 30 прогонов завершили 20 с; физические ошибки приведены без PASS-порога.

## Данные

| Файл | Назначение |
|---|---|
| [report.json](report.json) | 30 результатов, аудит MSE, команды, хэши, настройки и версии библиотек |
| [manifest.json](manifest.json) | Commit и исходные авторские файлы; статус MuJoCo; расположение трасс |
| [cadence.cpp](cadence.cpp), [cadence.json](cadence.json) | Потеря энергии при локальном дроблении даже без CCD |
| [bullet-small-step-repeat.json](bullet-small-step-repeat.json) | Повтор выброса Bullet при h=0.5 мс |

CSV и журналы оставлены в `/tmp/simbenchmark-final-bouncing-v2-20261001`;
SHA-256 каждого CSV есть в отчёте. Авторские исходники сохранены в
`/tmp/phys-ccd-retry-20261001`. Каталоги `/tmp` временны; для длительного
хранения их следует скопировать. Репозиторий содержит итоговые данные и
код воспроизведения без локальных бинарников и больших CSV.

## Воспроизведение на этой машине

Внешние библиотеки собраны отдельно; commits и хэши есть в отчёте.
SDK остаётся зависимым только от стандартной библиотеки и Threads.

```sh
cmake -S tools/rigid_reference -B /tmp/phys-engine-reference-build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="/tmp/phys-jolt-reference-install;/tmp/phys-box3d-reference-install"
cmake --build /tmp/phys-engine-reference-build --target rf_simbenchmark_bouncing --parallel 6
RF_THREADS=10 nice -n 10 taskset -c 1-5,7-11 \
  python3 tools/run_simbenchmark_bouncing.py \
  --binary /tmp/phys-engine-reference-build/rf_simbenchmark_bouncing \
  --upstream /tmp/phys-ccd-retry-20261001 \
  --provenance /tmp/simbenchmark-reference-provenance-20261001.json \
  --output /tmp/simbenchmark-NEW
python3 tools/plot_simbenchmark_bouncing.py /tmp/simbenchmark-NEW \
  docs/book/evidence/simbenchmark-bouncing-20261001
```

`NEW` должен быть новым каталогом. Процессы идут последовательно, с пределом
180 с на конфигурацию и проверкой неизменности бинарника. Для ошибки и
времени используются разные проходы. Отказ и тайм-аут сохраняются;
неполный опыт не получает оценки. `--engines sdk bullet` позволяет работать
без Jolt/Box3D. `--provenance` необязателен; на другой машине нужен свой
receipt с версиями и хэшами библиотек.

Срез Release, `RF_WERROR=ON`, `RF_STRICT_FP=OFF`; это не проверка строгого
FP-профиля CI. Контрольный опыт:

```sh
c++ -std=c++17 -O3 -Wall -Wextra -Werror -Isrc -I. -isystem /usr/include/bullet \
  docs/book/evidence/simbenchmark-bouncing-20261001/cadence.cpp \
  build-core/librfcore.a -lBulletDynamics -lBulletCollision -lLinearMath -pthread \
  -o /tmp/rf_bouncing_cadence
RF_THREADS=10 nice -n 10 taskset -c 1-5,7-11 /tmp/rf_bouncing_cadence
```

Контрольный опыт повторно использует наш адаптер. Наблюдатель и аналитическая
формула не вызывают helper энергии SDK; это не независимая внешняя
реализация всего решателя.

## Ограничения

Авторский Bullet wrapper заменён дискретным API. У Box3D иной бюджет и
модель контакта; два профиля представлены отдельно. Времена отражают разные
бюджеты и потоки, общего рейтинга скорости нет. MuJoCo не запускался.
Проникновение измерено в принятых кадрах, непрерывного сертификата нет.
Пробные каталоги сохранены в `/tmp`, включая запуск во время пересборки;
итоговый отчёт использует только проверенный неизменный бинарник.
