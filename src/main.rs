use std::env;
use std::fs;

use anyhow::Result;
use anyhow::anyhow;

use lovers::config::Config;
use lovers::io::read_model;
use lovers::messages;

fn main() -> Result<()> {
    messages::display_logo();
    messages::display_lovers();

    // Collect command line arguments
    let args: Vec<String> = env::args().skip(1).collect();

    // Ensure that the user has provided a path to a TOML file
    if args.len() != 1 {
        eprintln!("Usage: lovers <path_to_config>");
        return Ok(());
    }

    // Try to get the absolute path and handle errors properly
    let path = match fs::canonicalize(&args[0]) {
        Ok(p) => p, // Successfully resolved to an absolute path
        Err(e) => {
            return Err(e.into());
        }
    };

    // Ensure the path exists after resolving it
    if !path.exists() {
        eprintln!("Config file does not exist");
        return Ok(());
    }

    let config = Config::load(
        path.to_str()
            .ok_or_else(|| anyhow!("Error: The canonicalized path is not valid UTF-8."))?,
    )?;

    read_model(config.models[0].as_str())?;

    Ok(())
}
