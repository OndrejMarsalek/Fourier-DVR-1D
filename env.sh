ENV_BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PYTHONPATH=$ENV_BASE_DIR:$PYTHONPATH

unset ENV_BASE_DIR

