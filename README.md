# parallels
Это репозиторий с задачами по параллельному программированию

# task1
Чтобы сгенерировать файлы сборки используйте одну из следующих команд:

- `cmake -B build` - сборка по умолчанию (с `float`)

- `cmake -B build -D USE_DOUBLE=ON` - сборка с `double`

- `cmake -B build -D USE_DOUBLE=OFF` - сборка с `float`

Сборка и запуск:

- `make -C build`

- `./build/sin_array_sum`
