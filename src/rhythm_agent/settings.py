"""Application settings loaded from environment and project configuration."""

from pydantic import Field, SecretStr
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Runtime settings; scientific protocols remain versioned separately."""

    environment: str = "development"
    data_dir: str = "data"
    azure_openai_endpoint: str | None = None
    azure_openai_api_key: SecretStr | None = None
    azure_openai_deployment: str | None = None
    azure_openai_timeout_seconds: float = Field(default=60.0, gt=0)
    azure_openai_max_retries: int = Field(default=3, ge=0)

    model_config = SettingsConfigDict(
        env_prefix="RHYTHM_AGENT_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

'''pydantic_setting类负责直接从env文件里面加载配置,一般就只有这个时候会用到'''